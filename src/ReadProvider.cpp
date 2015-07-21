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
static uint m_entryCount;
static std::map<uloc, float> iTable; // fallback

uint const estimateSize = 100;
//uint const estimateSize = 10000000;
uint const estimateStepSize = 1000;
uint const estimateThreshold = 1000;
uint const maxReadLength = 1000;

float * maxHitTable;
int maxHitTableIndex;
int m_CurrentReadLength;

ReadProvider::ReadProvider() :
		parser1(0), parser2(0), peDelimiter(
		Config.GetString("pe_delimiter")[0]), isPaired(
		Config.GetInt("paired") > 0), skipMateCheck(
				Config.GetInt(SKIP_MATE_CHECK) == 1) {

}

int CollectResultsFallback(int const readLength) {
	float maxCurrent = 0;

	for (std::map<uloc, float>::iterator itr = iTable.begin();
			itr != iTable.end(); itr++) {
		maxCurrent = std::max(maxCurrent, itr->second);
	}

	static const int skip = (
	Config.Exists("kmer_skip") ? Config.GetInt("kmer_skip", 0, -1) : 0) + 1;
	//float max = (seq->seq.l - CS::prefixBasecount + 1) / skip;

	int max = ceil((readLength - CS::prefixBasecount + 1) / skip * 1.0);

	if (max > 1.0f && maxCurrent <= max) {
		maxHitTable[maxHitTableIndex++] = (maxCurrent / ((max))); // * 0.85f + 0.05f;
	}

	iTable.clear();
	return 0;
}

static void PrefixSearch(ulong prefix, uloc pos, ulong mutateFrom,
		ulong mutateTo, void* data) {

	RefEntry const * entries = m_RefProvider->GetRefEntry(prefix, m_entry); // Liefert eine liste aller Vorkommen dieses Praefixes in der Referenz
	RefEntry const * cur = entries;

	for (int i = 0; i < m_entryCount; i++) {
		//Get kmer-weight.
//		float weight = cur->weight;
		float weight = 1.0f;

		int const n = cur->refCount;

		if (cur->reverse) {
			for (int i = 0; i < n; ++i) {
				uloc curLoc = cur->getRealLocation(cur->ref[i])
						- (m_CurrentReadLength - (pos + CS::prefixBasecount));
				curLoc = GetBin(curLoc); // position offset
				if (iTable.count(curLoc) == 0) {
					iTable[curLoc] = weight;
				} else {
					iTable[curLoc] += weight;
				}
			}

		} else {
			for (int i = 0; i < n; ++i) {
				uloc curLoc = cur->getRealLocation(cur->ref[i]) - pos;
				curLoc = GetBin(curLoc); // position offset
				if (iTable.count(curLoc) == 0) {
					iTable[curLoc] = weight;
				} else {
					iTable[curLoc] += weight;
				}
			}
		}

		cur++;
	}
}

void PrefixMutateSearchEx(ulong prefix, uloc pos, ulong mutateFrom,
		ulong mutateTo, void* data, int mpos = 0);

void PrefixMutateSearch(ulong prefix, uloc pos, ulong mutateFrom,
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

void PrefixMutateSearchEx(ulong prefix, uloc pos, ulong mutateFrom,
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
	typedef void (*PrefixIterationFn)(ulong prefix, uloc pos, ulong mutateFrom,
			ulong mutateTo, void* data);
	PrefixIterationFn fnc = &PrefixSearch;

	bool const isPaired = Config.GetInt("paired") > 1;

	char const * const fileName1 =
	Config.Exists("qry1") ?
	Config.GetString("qry1") :
							Config.GetString("qry");
	char const * const fileName2 =
	Config.Exists("qry2") ? Config.GetString("qry2") : 0;

	bool m_EnableBS = false;
	m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

	if (m_EnableBS)
		fnc = &PrefixMutateSearch;

	Log.Verbose("Initializing ReadProvider");

	uint readCount = 0;
	Timer tmr;
	tmr.ST();
	size_t maxLen = 0;
	bool estimate = !(Config.Exists("skip_estimate")
			&& Config.GetInt("skip_estimate"));
	if (!Config.Exists("qry_max_len") || estimate) {
		//default value for estimation
		static int const qryMaxLen = 10000;
		parser1 = DetermineParser(fileName1, qryMaxLen);
		if (estimate) {
			Log.Message("Estimating parameter from data");

			maxHitTable = new float[estimateSize];
			maxHitTableIndex = 0;

			m_RefProvider = NGM.GetRefProvider(0);

			m_entryCount = m_RefProvider->GetRefEntryChainLength();
			m_entry = new RefEntry[ m_entryCount ];
		}

		int l = 0;

		size_t minLen = 9999999;
		size_t sumLen = 0;

		bool finish = false;

		MappedRead * read = new MappedRead(0, qryMaxLen);

		try {
			while ((l = parser1->parseRead(read)) >= 0 && !finish) {
				if (l > 0) {
					maxLen = std::max(maxLen, (size_t) read->length);
					minLen = std::min(minLen, (size_t) read->length);
					sumLen += read->length;

//					Log.Message("Name: %s", read->name);
//					Log.Message("Read: %s", read->Seq);
//					Log.Message("Qlty: %s", read->qlty);

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
						CS::PrefixIteration((char const *) read->Seq, read->length, fnc, mutateFrom, mutateTo, (void *) this, (uint) 0, 0);
						CollectResultsFallback(m_CurrentReadLength);
					} else if (readCount == (estimateSize + 1)) {
						if ((maxLen - minLen) < 10 && !Config.Exists(ARGOS)) {
							finish = true;
						} else {
							if (Config.GetInt(RLENGTH_CHECK)) {
								Log.Warning("Input reads don't have the same length!");
								Log.Warning("Parameter 'force-rlength-check' found. Determining max. read length now. This might take some time!");
							} else {
								Log.Warning("Input reads don't have the same length!");
								Log.Warning("Maximum read length found in the first %d reads is %d. For longer reads only the first %d bp will be mapped.", estimateSize, maxLen, (int) (maxLen * 1.1f));
								maxLen = maxLen * 1.1f;
								Log.Warning("The maximum read length can be overwritten with the '--max-read-length' parameter. With '--force-rlength-check', NextGenMap will run through all reads to find the max. read length. This might take some time.");
								finish = true;
							}
						}
					}
				}
			}
		} catch(char * ex) {
			Log.Error("%s", ex);
			Fatal();
		}
		delete read;
		read = 0;
		if (!finish) {
			Log.Message("Reads found in files: %d", readCount);
			NGM.Stats->TotalSeqs = readCount;
		}
		if (readCount == 0) {
			Log.Error("No reads found in input file.");
			Fatal();
		}
//		if (Config.GetInt(MAX_READ_LENGTH) > 0) {
//			Log.Warning("Max. read length overwritten bei user: %d", Config.GetInt(MAX_READ_LENGTH));
//			maxLen = Config.GetInt(MAX_READ_LENGTH);
//		}
		maxLen = (maxLen | 1) + 1;
		int avgLen = sumLen / readCount;

//		if (maxLen > maxReadLength) {
//			Log.Warning("Max. supported read length is 1000bp. All reads longer than 1000bp will be hard clipped in the output.");
//			maxLen = maxReadLength;
//		}

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
			for (int i = 0; i < maxHitTableIndex; ++i) {
				sum += maxHitTable[i];
			}

			float m_CsSensitivity = 0.0f;
			if (!m_EnableBS) {
				static const int skip = (
				Config.Exists("kmer_skip") ?
				Config.GetInt("kmer_skip", 0, -1) :
												0) + 1;
				//float max = (avgLen - CS::prefixBasecount + 1) / skip;
				int max = ceil((avgLen - CS::prefixBasecount + 1) / skip * 1.0);
				float avg = sum / maxHitTableIndex * 1.0f;

				float avgHit = max * avg;

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
			}
			delete[] maxHitTable;
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
	Log.Message("Estimating parameter took %.3fs", tmr.ET());

	parser1 = DetermineParser(fileName1, maxLen);
	if (fileName2 != 0) {
		parser2 = DetermineParser(fileName2, maxLen);
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

	int l = 0;

	try {
		l = parser->parseRead(read);

		if (l >= 0) {
			int nameLength = strlen(read->name);

			if (isPaired && read->name[nameLength - 2] == peDelimiter) {
				nameLength -= 2;
				read->name[nameLength] = '\0';
			}

			Log.Debug(2, "READ_%d\tINPUT\t%s", id, read->name);
			Log.Debug(16384, "READ_%d\tINPUT_DETAILS\t%s\t%s\t%s\t%s", id, read->Seq, read->qlty, read->AdditionalInfo);

			NGM.AddReadRead(read->ReadId);
		} else {

			Log.Debug(2, "READ_%d\tINPUT\t%s error while reading", id, read->name);

			if(l == -2) {
				Log.Error("Read %s: Length of read not equal length of quality values.", read->name);
				Fatal();
			} else if (l != -1) {
				//TODO correct number when paired
				Log.Error("Unknown error while parsing read number %d (error code: %d)", id + 1, l);
				Fatal();
			}
			delete read;
			read = 0;
		}
	} catch (char * ex) {
		Log.Error("%s", ex);
		Fatal();
	}
	return read;
}

IParser * ReadProvider::DetermineParser(char const * fileName,
		int const qryMaxLen) {

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
		parser = new SamParser(qryMaxLen);
	} else {
		if (strncmp(buffer, "BAM", 3) == 0) {
#ifdef _BAM
			Log.Message("Input is BAM");
			parser= new BamParser(qryMaxLen);
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
			parser = new FastXParser(qryMaxLen);
		}
	}
	gzclose(fp);
	delete[] buffer;
	parser->init(fileName, Config.GetInt(KEEPTAGS) == 1);
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
			read2 = NextRead(parser2, readid2);
		}

		if (read1 != 0 && read2 != 0) {
			if (!skipMateCheck && strcmp(read1->name, read2->name) != 0) {
				Log.Error("Error while reading paired end reads.");
				Log.Error("Names of mates don't match: %s and %s.", read1->name, read2->name);
				Log.Error("NextGenMap expects paired end read names with the format: <read name>/<mate nunber> (e.g. @HWI-ST1176_0172:8:1101:1234:1934/1)");
				Log.Error("Use -d/--pe-delimiter to specify a different delimiter than '/'");
				Log.Error("If the format of the read names is correct this error might be caused by missing mates in the input file. Please check your input files or use ngm-utils interleave to match mate pairs.");
				Log.Error("If you are sure that your input files are valid use --skip-mate-check");
				Fatal();
			}
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
			delete read->Paired;
			delete read;
		} else {
			// Paired read is still in use...mark read for deletion
			read->SetFlag(NGMNames::DeletionPending);
		}
	} else {
		// Single mode or no existing pair
		delete read;
	}
}
