/*
 * ReadProvider.cpp
 *
 *  Created on: Jun 14, 2012
 *      Author: philipp_
 */

#include "ReadProvider.h"

#include <zlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <limits.h>

#include "Config.h"
#include "Log.h"
#include "Timing.h"
#include "CS.h"
#include "MappedRead.h"
#include "FastxParser.h"
#include "SamParser.h"

#ifdef _BAM
#include "BamParser.h"
#endif

using NGMNames::ReadStatus;

#undef module_name
#define module_name "INPUT"

static IRefProvider const * m_RefProvider;
static RefEntry * m_entry;
static std::map<SequenceLocation, float> iTable; // fallback

//uint const estimateSize = 10000;
uint const estimateSize = 10000000;
uint const estimateStepSize = 1000;
uint const estimateThreshold = 1000;
uint const maxReadLength = 1000;

float * maxHitTable;
#ifdef _DEBUGRP
float * maxHitTableDebug;
#endif
int maxHitTableIndex;
int m_CurrentReadLength;

ReadProvider::ReadProvider() :
		parser1(0), parser2(0), peDelimiter(
				Config.GetString("pe_delimiter")[0]), isPaired(
				Config.GetInt("paired") > 0) {

}

inline int GetBin(uint pos) {
	static int shift = CS::calc_binshift(20);
	return pos >> shift;
	//return pos;
}

int CollectResultsFallback(int const readLength) {
	float maxCurrent = 0;

	Log.Verbose("-maxHitTableIndex: %d", maxHitTableIndex);
	for (std::map<SequenceLocation, float>::iterator itr = iTable.begin();
			itr != iTable.end(); itr++) {
		maxCurrent = std::max(maxCurrent, itr->second);

		//if(maxHitTableIndex == 8) {
		Log.Verbose("---Key: %d %u %d, maxCurrent: %f, itr->second: %f", itr->first.getrefId(), itr->first.m_Location, (int)itr->first.isReverse(), maxCurrent, itr->second);
		//}
	}

	static const int skip = (
			Config.Exists("kmer_skip") ? Config.GetInt("kmer_skip", 0, -1) : 0)
			+ 1;
	//float max = (seq->seq.l - CS::prefixBasecount + 1) / skip;

	int max = ceil((readLength - CS::prefixBasecount + 1) / skip * 1.0);

	if (max > 1.0f && maxCurrent <= max) {
#ifdef _DEBUGRP
		maxHitTableDebug[maxHitTableIndex] = maxCurrent;
#endif
		maxHitTable[maxHitTableIndex++] = (maxCurrent / ((max))); // * 0.85f + 0.05f;
		//Log.Message("Result: %f, %f, %f, %f -> %f", maxCurrent, maxCurrent / seq->seq.l, max, max / seq->seq.l, maxHitTable[maxHitTableIndex-1]);
	}

	iTable.clear();

	return 0;
}

static void PrefixSearch(ulong prefix, uint pos, ulong mutateFrom,
		ulong mutateTo, void* data) {

	RefEntry const * cur = m_RefProvider->GetRefEntry(prefix, m_entry); // Liefert eine liste aller Vorkommen dieses Praefixes in der Referenz

	while (cur != 0) {
		//Get kmer-weight.
//		float weight = cur->weight;
		float weight = 1.0f;

		int const n = cur->refCount;

		Log.Verbose("------Prefix: %u, RefCount: %d", prefix, n);

		if (cur->reverse) {
			for (int i = 0; i < n; ++i) {
				SequenceLocation curLoc = cur->ref[i];
				curLoc.setRefId(1);
				curLoc.setReverse(true);
				curLoc.m_Location = GetBin(
						curLoc.m_Location
								- (m_CurrentReadLength
										- (pos + CS::prefixBasecount))); // position offset
				iTable[curLoc] += weight;
			}

		} else {
			for (int i = 0; i < n; ++i) {
				SequenceLocation curLoc = cur->ref[i];
				curLoc.setRefId(0);
				curLoc.setReverse(false);
				curLoc.m_Location = GetBin(curLoc.m_Location - pos); // position offset
				iTable[curLoc] += weight;
			}
		}

		cur = cur->nextEntry;
	}
}

void PrefixMutateSearchEx(ulong prefix, uint pos, ulong mutateFrom,
		ulong mutateTo, void* data, int mpos = 0);

void PrefixMutateSearch(ulong prefix, uint pos, ulong mutateFrom,
		ulong mutateTo, void* data) {
	static int const cMutationLocLimit =
			Config.Exists("bs_cutoff") ? Config.GetInt("bs_cutoff") : 6;
	ulong const mask = 0x3;

	int mutationLocs = 0;
	for (int i = 0; i < (int) CS::prefixBasecount; ++i) {
		ulong base = mask & (prefix >> (i * 2));
		if (base == mutateFrom)
			++mutationLocs;
	}

	if (mutationLocs <= cMutationLocLimit)
		PrefixMutateSearchEx(prefix, pos, mutateFrom, mutateTo, data);
}

void PrefixMutateSearchEx(ulong prefix, uint pos, ulong mutateFrom,
		ulong mutateTo, void* data, int mpos) {
	PrefixSearch(prefix, pos, mutateFrom, mutateTo, data);

	ulong const mask = 0x3;
	for (int i = mpos; i < (int) CS::prefixBasecount; ++i) {
		ulong cur = mask & (prefix >> (i * 2));

		if (cur == mutateFrom) {
			ulong p1 = (prefix & ~(mask << (i * 2)));
			ulong p2 = (mutateTo << (i * 2));
			PrefixMutateSearchEx(p1 | p2, pos, mutateFrom, mutateTo, data,
					i + 1);
		}
	}
}

uint ReadProvider::init() {
	typedef void (*PrefixIterationFn)(ulong prefix, uint pos, ulong mutateFrom,
			ulong mutateTo, void* data);
	PrefixIterationFn fnc = &PrefixSearch;

	bool const isPaired = Config.GetInt("paired") > 1;

	char const * const fileName1 =
			Config.Exists("qry1") ?
					Config.GetString("qry1") : Config.GetString("qry");
	char const * const fileName2 =
			Config.Exists("qry2") ? Config.GetString("qry2") : 0;

	bool m_EnableBS = false;
	m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

	if (m_EnableBS)
		fnc = &PrefixMutateSearch;

	Log.Message("Initializing ReadProvider");

	uint readCount = 0;
	Timer tmr;
	tmr.ST();
	size_t maxLen = 0;
	bool estimate = !(Config.Exists("skip_estimate")
			&& Config.GetInt("skip_estimate"));
	if (!Config.Exists("qry_max_len") || estimate) {
		parser1 = DetermineParser(fileName1);
		if (estimate) {
			Log.Message("Estimating parameter from data");

#ifdef _DEBUGRP
			maxHitTableDebug = new float[estimateSize];
#endif
			maxHitTable = new float[estimateSize];
			maxHitTableIndex = 0;

			m_RefProvider = NGM.GetRefProvider(0);
			m_entry = new RefEntry(0);
			m_entry->nextEntry = new RefEntry(0);
		}

		int l = 0;

		size_t minLen = 9999999;
		size_t sumLen = 0;

		bool finish = false;

		int const qryMaxLen = 10000;

		MappedRead * read = new MappedRead(0, qryMaxLen);

		while ((l = parser1->parseRead(read)) >= 0 && !finish) {
			if (l > 0) {
				maxLen = std::max(maxLen, (size_t) read->length);
				minLen = std::min(minLen, (size_t) read->length);
				sumLen += read->length;

//				Log.Message("Name: %s", read->name);
//				Log.Message("Read: %s", read->Seq);
//				Log.Message("Qlty: %s", read->qlty);

				readCount += 1;
				if (estimate && (readCount % estimateStepSize) == 0
						&& readCount < estimateSize) {
					ulong mutateFrom;
					ulong mutateTo;
					if (isPaired && (readCount & 1)) {
						//Second mate
						mutateFrom = 0x0;
						mutateTo = 0x3;
					} else {
						//First mate
						mutateFrom = 0x2;
						mutateTo = 0x1;
					}
					m_CurrentReadLength = read->length;
					Log.Verbose("-Iteration");
					CS::PrefixIteration((char const *) read->Seq,
							(uint) read->length, fnc, mutateFrom, mutateTo,
							(void *) this, (uint) 0, (uint) 0);
					Log.Verbose("-Collect: %s", parser1->read->name.s);
					CollectResultsFallback(m_CurrentReadLength);
				} else if (readCount == (estimateSize + 1)) {
					if ((maxLen - minLen) < 10) {
						finish = true;
					} else {
						Log.Warning("Reads don't have the same length. Determining max. read length now.");
					}
				}
			}
		}
		delete read;
		read = 0;
		if (!finish) {
			Log.Message("Reads found in files: %d", readCount);
		}
		if (readCount == 0) {
			Log.Error("No reads found in input file.");
			Fatal();
		}
		maxLen = (maxLen | 1) + 1;
		int avgLen = sumLen / readCount;

		if (maxLen > maxReadLength) {
			Log.Warning("Max. supported read length is 1000bp. All reads longer than 1000bp will be hard clipped in the output.");
			maxLen = maxReadLength;
		}
		((_Config*) _config)->Override("qry_max_len", (int) maxLen);
		((_Config*) _config)->Override("qry_avg_len", (int) avgLen);

		Log.Message("Average read length: %d (min: %d, max: %d)", avgLen, minLen, maxLen);

		if (Config.Exists("corridor")) {
			Log.Warning("Corridor witdh overwritten!");
		} else {
			((_Config*) _config)->Default("corridor", (int) (5 + avgLen * 0.15));
		}
		Log.Message("Corridor width: %d", Config.GetInt("corridor"));

		delete parser1;

		if (estimate && readCount >= estimateThreshold) {
			float sum = 0.0f;
#ifdef _DEBUGRP
			FILE* ofp;
			ofp = fopen((std::string(Config.GetString("output")) + std::string(".b")).c_str(), "w");
#endif
			for (int i = 0; i < maxHitTableIndex; ++i) {
				Log.Verbose("%f += %f", sum, maxHitTable[i]);
				sum += maxHitTable[i];
#ifdef _DEBUGRP
				fprintf(ofp, "%f\n", maxHitTableDebug[i]);
#endif
			}
#ifdef _DEBUGRP
			fclose(ofp);
#endif

			float m_CsSensitivity = 0.0f;
			if (!m_EnableBS) {
				static const int skip = (
						Config.Exists("kmer_skip") ?
								Config.GetInt("kmer_skip", 0, -1) : 0) + 1;
				//float max = (avgLen - CS::prefixBasecount + 1) / skip;
				int max = ceil((avgLen - CS::prefixBasecount + 1) / skip * 1.0);
				float avg = sum / maxHitTableIndex * 1.0f;

				float avgHit = max * avg;
				Log.Verbose("max: %d, sum: %f, maxHitTableIndex: %d, avg: %f, avgHit: %f", max, sum, maxHitTableIndex, avg, avgHit);

				Log.Message("Average kmer hits pro read: %f", avgHit);
				Log.Message("Max possible kmer hit: %d", max);
				//Log.Message("Max possible kmer hit (old): %f", (avgLen - CS::prefixBasecount + 1) / skip);

				//m_CsSensitivity = (avg / ((max / avgLen))) * 0.90f + 0.05f;
				m_CsSensitivity = std::min(std::max(0.3f, avg), 0.9f);
				Log.Green("Estimated sensitivity: %f", m_CsSensitivity);

				if (Config.Exists("sensitivity")) {
					//float x = Config.GetFloat("sensitivity", 0, 100);
					//m_CsSensitivity = (100 - x) / 100.0f;
					m_CsSensitivity = Config.GetFloat("sensitivity", 0, 1);
					Log.Warning("Sensitivity threshold overwritten by user. Using %f", m_CsSensitivity);
				}
				((_Config*) _config)->Override("sensitivity", m_CsSensitivity);
			} /*else {
			 m_CsSensitivity = 0.5f;
			 if (Config.Exists("sensitivity")) {
			 m_CsSensitivity = Config.GetFloat("sensitivity", 0, 1);
			 }
			 }*/

			//std::sort(lengthTable, lengthTable + lengthTableIndex);
			//int size = lengthTable[(int) (0.5f * lengthTableIndex) - 1] * 32;
			//int bits = (int) std::min(20.0, ceil(log2(size)));
			//Log.Green("Optimale hash table size: %d (%d)", (int)pow(2.0, (double)bits), lengthTable[lengthTableIndex - 1]);
			//if(Config.Exists("cs_tablen")) {
			//	bits = Config.GetInt("cs_tablen", 0, 32);
			//	Log.Warning("SearchTableLength overwritten by user. Using %d bits", (int)pow(2.0, (double)bits));
			//}
			//((_Config*) _config)->Override("cs_tablen", bits);
			delete[] maxHitTable;
#ifdef _DEBUGRP
			delete[] maxHitTableDebug;
#endif
		} else {
			Log.Warning("Not enough reads to estimate parameter");
		}

	} else {
		Log.Warning("qry_max_len was found in config file. Please make sure that this number is correct: %d", Config.GetInt("qry_max_len"));
	}

	if (!Config.Exists("sensitivity")) {
		((_Config*) _config)->Override("sensitivity", 0.5f);
		if (!m_EnableBS) {
			Log.Warning("Sensitivity parameter neither set nor estimated. Falling back to default.");
		} else {
			Log.Message("Sensitivity parameter set to 0.5");
		}
	}
	Log.Message("Initializing took %.3fs", tmr.ET());

	parser1 = DetermineParser(fileName1);
	if (fileName2 != 0) {
		parser2 = DetermineParser(fileName2);
	}

	return 0;
}

ReadProvider::~ReadProvider() {
	if (parser1 != 0) {
		delete parser1;
		parser1 = 0;
	}
	if (parser2 != 0) {
		delete parser2;
		parser2 = 0;
	}
}

MappedRead * ReadProvider::NextRead(IParser * parser, int const id) {
	static int const qryMaxLen = Config.GetInt("qry_max_len");
	MappedRead * read = new MappedRead(id, qryMaxLen);

	int l = parser->parseRead(read);

	if (l > 0) {
		int nameLength = strlen(read->name);

		if (isPaired && read->name[nameLength - 2] == peDelimiter) {
			nameLength -= 2;
			read->name[nameLength] = '\0';
		}

//		Log.Message("Name: %s", read->name);
//		Log.Message("Read: %s", read->Seq);
//		Log.Message("Qlty: %s", read->qlty);
//		if(read->AdditionalInfo != 0) {
//			Log.Message("Ainf: %s", read->AdditionalInfo);
//		} else {
//			Log.Message("No additional info.");
//		}

		NGM.AddReadRead(read->ReadId);
	} else {
		if(l == -2) {
			Log.Error("Read %s: Length of read not equal length of quality values.", read->name);
			Fatal();
		} else if (l != -1) {
			//TODO correct number when paired
			Log.Error("Unknown error while parsing read %d (%d)", id + 1, l);
			if(isPaired) {
				Fatal();
			}
		}
		delete read;
		read = 0;
	}
//	if (l >= 0) {
//		if (parser->read->seq.l == parser->read->qual.l || parser->read->qual.l == 0) {
//			read = new MappedRead(id, qryMaxLen);
//
//			//Name
//			static size_t const MAX_READNAME_LENGTH = 100;
//			read->name = new char[MAX_READNAME_LENGTH];
//			int nameLength = std::min(MAX_READNAME_LENGTH - 1, parser->read->name.l);
//
//			if(isPaired && parser->read->name.s[nameLength - 2] == peDelimiter) {
//				nameLength -= 2;
//			}
//
//			memcpy(read->name, parser->read->name.s, nameLength);
//			read->name[nameLength] = '\0';
//
//			//Sequence
//			read->Seq = new char[qryMaxLen];
//			memset(read->Seq, '\0', qryMaxLen);
//			if (parser->read->seq.l != 0) {
//				read->length = std::min(parser->read->seq.l, (size_t) qryMaxLen - 1);
//				int nCount = 0;
//				for (int i = 0; i < read->length; ++i) {
//					char c = toupper(parser->read->seq.s[i]);
//					if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
//						read->Seq[i] = c;
//					} else {
//						read->Seq[i] = 'N';
//						nCount += 1;
//					}
//
//				}
//			} else {
//				Log.Verbose("Empty read found (%s). Filling with Ns.", read->name);
//				read->length = qryMaxLen - 2;
//				memset(read->Seq, 'N', read->length);
//				read->SetFlag(NGMNames::Empty);
//			}
//
//			//Quality
//			read->qlty = 0;
//			if (parser->read->qual.l > 0) {
//				read->qlty = new char[read->length + 1];
//				memcpy(read->qlty, parser->read->qual.s, read->length);
//				read->qlty[read->length] = '\0';
//			}
//
//			NGM.AddReadRead(read->ReadId);
//		} else {
//			Log.Error("Discarding read %s. Length of read not equal length of quality values.", parser->read->name.s);
//			Fatal();
//		}
//	} else {
//		if (l == -1) {
//			Log.Verbose("End of input file reached.");
//		} else {
//			Log.Error("Error while parsing read %d (%d)", id, l);
//			Fatal();
//		}
//	}
	return read;
}

IParser * ReadProvider::DetermineParser(char const * fileName) {
	gzFile fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("File %s does not exist!",fileName);
		Fatal();
	}
	IParser * parser = 0;
	char * buffer = new char[1000];
	while (gzgets(fp, buffer, 1000) > 0 && buffer[0] == '@') {
	}

	int count = 0;
	for (size_t i = 0; i < 1000 && buffer[i] != '\0' && buffer[i] != '\n';
			i++) {
		if (buffer[i] == '\t') {
			count++;
		}
	}
	if (count >= 10) {
		Log.Message("Input is SAM");
		parser = new SamParser();
	} else {
		if (strncmp(buffer, "BAM", 3) == 0) {
#ifdef _BAM
			Log.Message("Input is BAM");
			parser= new BamParser();
#else
			Log.Error("BAM input detected. NGM was compiled without BAM support!");
			Fatal();
#endif
		} else {
			if (buffer[0] == '>') {
				Log.Message("Input is Fasta");
			} else {
				Log.Message("Input is Fastq");
			}
			parser = new FastXParser();
		}
	}
	gzclose(fp);
	delete[] buffer;
	parser->init(fileName);
	return parser;
}

MappedRead * ReadProvider::GenerateSingleRead(int const readid) {
	MappedRead * read = 0;

	read = NextRead(parser1, readid); //

	return read;
}

// Sequential (important for pairs!) read generation
bool ReadProvider::GenerateRead(int const readid1, MappedRead * & read1,
		int const readid2, MappedRead * & read2) {

	if (isPaired) {
		static bool const isInterleaved = Config.Exists("qry");
		if (isInterleaved) {
			read1 = GenerateSingleRead(readid1);
			read2 = GenerateSingleRead(readid2);
		} else {
			read1 = NextRead(parser1, readid1);
			read2 = NextRead(parser2, readid1);
		}

		if (read1 != 0 && read2 != 0) {
			read1->Paired = read2;
			read2->Paired = read1;
			return true;
		} else if (read1 == 0 && read2 == 0) {
			return false;
		} else {
			Log.Error("Error in input file. Number of reads in input not even. Please check the input or mapped in single-end mode.");
			Fatal();
		}
	} else {
		read1 = GenerateSingleRead(readid1);
		read2 = GenerateSingleRead(readid2);
		return read1 != 0 && read2 != 0;
	}
}

void ReadProvider::DisposeRead(MappedRead * read) {
	static bool const isPaired = Config.GetInt("paired") > 0;
	if (isPaired) // Program runs in paired mode and pair exists
	{
		if (read->Paired->HasFlag(NGMNames::DeletionPending)) // Paired read was marked for deletion
				{
			// delete paired and read
			Log.Verbose("Deleting read %d (%s)", read->ReadId, read->name);
			Log.Verbose("Deleting read %d (%s)", read->Paired->ReadId, read->Paired->name);
			delete read->Paired;
			delete read;
		} else {
			// Paired read is still in use...mark read for deletion
			read->SetFlag(NGMNames::DeletionPending);
		}
	} else {
		// Single mode or no existing pair
		Log.Verbose("Deleting single read %d (%s)", read->ReadId, read->name);
		delete read;
	}
}
