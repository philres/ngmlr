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

#include "kseq.h"
#include "Config.h"
#include "Log.h"
#include "Timing.h"
#include "CS.h"
#include "MappedRead.h"
#include "IParser.h"
#include "FastxParser.h"
#include "SamParser.h"

#ifdef _BAM
#include "BamParser.h"
#endif

using NGMNames::ReadStatus;

#undef module_name
#define module_name "READPROV"

kseq_t *seq;

IParser * parser;

static IRefProvider const * m_RefProvider;
static RefEntry * m_entry;
static std::map<SequenceLocation, float> iTable; // fallback

//uint const estimateSize = 10000;
uint const estimateSize = 1000000;
uint const estimateStepSize = 1000;
uint const estimateThreshold = 1000;
uint const maxReadLength = 1000;

float * maxHitTable;
int maxHitTableIndex;
int m_CurrentReadLength;

ReadProvider::ReadProvider() {

}

inline int GetBin(uint pos) {
	static int shift = CS::calc_binshift(20);
	return pos >> shift;
	//return pos;
}

int CollectResultsFallback() {
	float maxCurrent = 0;

	for (std::map<SequenceLocation, float>::iterator itr = iTable.begin(); itr != iTable.end(); itr++) {
		maxCurrent = std::max(maxCurrent, itr->second);
	}

	static const int skip = (Config.Exists("kmer_skip") ? Config.GetInt("kmer_skip", 0, -1) : 0) + 1;
	//float max = (seq->seq.l - CS::prefixBasecount + 1) / skip;

	int max = ceil((seq->seq.l - CS::prefixBasecount + 1) / skip * 1.0);


	if (max > 1.0f && maxCurrent <= max ) {
		maxHitTable[maxHitTableIndex++] = (maxCurrent / ((max)));// * 0.85f + 0.05f;
		//Log.Message("Result: %f, %f, %f, %f -> %f", maxCurrent, maxCurrent / seq->seq.l, max, max / seq->seq.l, maxHitTable[maxHitTableIndex-1]);
	}

	iTable.clear();

	return 0;
}

static void PrefixSearch(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data) {

	RefEntry const * cur = m_RefProvider->GetRefEntry(prefix, m_entry); // Liefert eine liste aller Vorkommen dieses Praefixes in der Referenz

	while (cur != 0) {
		//Get kmer-weight.
//		float weight = cur->weight;
		float weight = 1.0f;

		int const n = cur->refCount;

		if (cur->reverse) {
			for (int i = 0; i < n; ++i) {
				SequenceLocation curLoc = cur->ref[i];
				curLoc.m_RefId = 1;
				curLoc.m_Location = GetBin(curLoc.m_Location - (m_CurrentReadLength - (pos + CS::prefixBasecount))); // position offset
				iTable[curLoc] += weight;
			}

		} else {
			for (int i = 0; i < n; ++i) {
				SequenceLocation curLoc = cur->ref[i];
				curLoc.m_RefId = 0;
				curLoc.m_Location = GetBin(curLoc.m_Location - pos); // position offset
				iTable[curLoc] += weight;
			}
		}

		cur = cur->nextEntry;
	}
}

void PrefixMutateSearchEx(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data, int mpos = 0);

void PrefixMutateSearch(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data) {
	static int const cMutationLocLimit = Config.Exists("bs_cutoff") ? Config.GetInt("bs_cutoff") : 6;
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

void PrefixMutateSearchEx(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data, int mpos) {
	PrefixSearch(prefix, pos, mutateFrom, mutateTo, data);

	ulong const mask = 0x3;
	for (int i = mpos; i < (int) CS::prefixBasecount; ++i) {
		ulong cur = mask & (prefix >> (i * 2));

		if (cur == mutateFrom) {
			ulong p1 = (prefix & ~(mask << (i * 2)));
			ulong p2 = (mutateTo << (i * 2));
			PrefixMutateSearchEx(p1 | p2, pos, mutateFrom, mutateTo, data, i + 1);
		}
	}
}

uint ReadProvider::init(char const * fileName) {
	typedef void (*PrefixIterationFn)(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	PrefixIterationFn fnc = &PrefixSearch;

	bool m_EnableBS = false;
//	if (Config.Exists("bs_mapping"))
	m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

	if (m_EnableBS)
		fnc = &PrefixMutateSearch;

	Log.Message("Initializing ReadProvider");

	uint readCount = 0;
	Timer tmr;
	tmr.ST();
	size_t maxLen = 0;
	bool estimate = !(Config.Exists("skip_estimate") && Config.GetInt("skip_estimate"));
	if (!Config.Exists("qry_max_len") || estimate) {
		DetermineParser(fileName);
		if (estimate) {
			Log.Message("Estimating parameter from data");

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
		while ((l = parser->parseRead(seq)) > 0 && !finish) {
			maxLen = std::max(maxLen, seq->seq.l);
			minLen = std::min(minLen, seq->seq.l);
			sumLen += seq->seq.l;

			//Log.Message("Name: %s", seq->name.s);
			//Log.Message("Read: %s", seq->seq.s);

			readCount += 1;
			if (estimate && (readCount % estimateStepSize) == 0 && readCount < estimateSize) {
				ulong mutateFrom;
				ulong mutateTo;
				if (NGM.Paired() && (readCount & 1)) {
					//Second mate
					mutateFrom = 0x0;
					mutateTo = 0x3;
				} else {
					//First mate
					mutateFrom = 0x2;
					mutateTo = 0x1;
				}
				m_CurrentReadLength = seq->seq.l;
				CS::PrefixIteration((char const *) seq->seq.s, (uint) seq->seq.l, fnc, mutateFrom, mutateTo, (void *) this, (uint) 0,
						(uint) 0);
				CollectResultsFallback();
			} else if (readCount == (estimateSize + 1)) {
				if ((maxLen - minLen) < 10) {
					finish = true;
				} else {
					Log.Warning("454/Ion Torrent data set detected. Determining max. read length now.");
				}
			}
		}
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
		((_Config*) _config)->Default("corridor", (int) (5 + avgLen * 0.15));
		Log.Message("Average read length: %d (min: %d, max: %d)", avgLen, minLen, maxLen);
		Log.Message("Corridor width: %d", Config.GetInt("corridor"));

		kseq_destroy(seq);
		delete parser;

		if (estimate && readCount >= estimateThreshold) {
			float sum = 0.0f;
			for (int i = 0; i < maxHitTableIndex; ++i) {
				sum += maxHitTable[i];
			}

			float m_CsSensitivity = 0.0f;
			if (!m_EnableBS) {
				static const int skip = (Config.Exists("kmer_skip") ? Config.GetInt("kmer_skip", 0, -1) : 0) + 1;
				//float max = (avgLen - CS::prefixBasecount + 1) / skip;
				int max = ceil((avgLen - CS::prefixBasecount + 1) / skip * 1.0);
				float avg = sum / maxHitTableIndex * 1.0f;

				float avgHit = max * avg;

				Log.Message("Average kmer hits pro read: %f", avgHit);
				Log.Message("Max possible kmer hit: %d", max);
				Log.Message("Max possible kmer hit (old): %f", (avgLen - CS::prefixBasecount + 1) / skip);

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
	DetermineParser(fileName);


	return 0;
}

ReadProvider::~ReadProvider() {
	delete parser;
	kseq_destroy(seq);
}

MappedRead * ReadProvider::NextRead(int const id) {
//	Log.Message("next REad");
	static int const qryMaxLen = Config.GetInt("qry_max_len");
	MappedRead * read = 0;
	int l = parser->parseRead(seq);
//	Log.Message("Name (%d): %s", seq->name.l, seq->name.s);
//	Log.Message("Seq  (%d): %s", seq->seq.l, seq->seq.s);
//	Log.Message("Qual (%d): %s", seq->qual.l, seq->qual.s);
	if (l >= 0) {
		if (seq->seq.l == seq->qual.l || seq->qual.l == 0) {
			read = new MappedRead(id);

			//Name
			static size_t const MAX_READNAME_LENGTH = 100;
			read->name = new char[MAX_READNAME_LENGTH];
			read->nameLength = std::min(MAX_READNAME_LENGTH - 1, seq->name.l);
			memcpy(read->name, seq->name.s, read->nameLength);
			read->name[read->nameLength] = '\0';

//			char const * debugRead = "FCC01PDACXX:4:1101:10342:37018#0/1";
//			if(strcmp(read->name, debugRead) == 0) {
//				Log.Error("Read %s found: assigning id %d", debugRead, read->ReadId);
//			}

			//Sequence
			read->Seq = new char[qryMaxLen];
			memset(read->Seq, '\0', qryMaxLen);
			read->length = std::min(seq->seq.l, (size_t) qryMaxLen - 1);
			int nCount = 0;
			for (int i = 0; i < read->length; ++i) {
				char c = toupper(seq->seq.s[i]);
				if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
					read->Seq[i] = c;
				} else {
					read->Seq[i] = 'N';
					nCount += 1;
				}

			}
//			if (nCount > qryMaxLen * 0.5f) {
//				Log.Warning("Discarding read %s (too many Ns)", read->name);
//				delete read;
//				return 0;
//			}
//			for (int i = read->length; i < qryMaxLen; ++i) {
//				read->Seq[i] = '\0';
//			}

			//memcpy(read->Seq, seq->seq.s, read->length);

			//Quality
			read->qlty = 0;
			if (seq->qual.l > 0) {
				read->qlty = new char[read->length + 1];
				memcpy(read->qlty, seq->qual.s, read->length);
				read->qlty[read->length] = '\0';
			}

//			Log.Message("%s", read->name);
//			Log.Message("%s", read->Seq);
//			if (read->qlty != 0)
//				Log.Message("%s", read->qlty);
			NGM.AddReadRead(read->ReadId);
		} else {
			Log.Error("Discarding read %s. Length of read not equal length of quality values.", seq->name.s);
			Fatal();
		}
	} else {
		if (l == -1) {
			Log.Message("End of input file reached.");
		} else {
			Log.Error("Error while parsing read %d (%d)", id, l);
			Fatal();
		}
	}
	return read;
}

//bool _SequenceProvider::IsValid(int n) const
//{
//	CheckQryNr(n);
//
//	if (NGM.Paired())
//	{
//		bool offstrand = (n & 1);
//		n >>= 1;
//
//		if (offstrand)
//			return prdIdx[n].Flags & 0x1;
//	}
//
//	return qryIdx[n].Flags & 0x1;
//}
void ReadProvider::DetermineParser(char const * fileName) {
	gzFile fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("File does not exist ",fileName);
	}
	char * buffer = new char[1000];
	while (gzgets(fp, buffer, 1000) > 0 && buffer[0] == '@') {
	}

	int count = 0;
	for (size_t i = 0; i < 1000 && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
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
	delete[] buffer;
	seq = parser->init_seq(fileName);
}

MappedRead * ReadProvider::GenerateSingleRead(int const readid) {
	MappedRead * read = 0;

//	if (SequenceProvider.IsValid(readid))
	read = NextRead(readid); //

	return read;
}

// Sequential (important for pairs!) read generation
MappedRead * ReadProvider::GenerateRead(int const readid) {
	static MappedRead * nextRead = 0;

	if (NGM.Paired()) {
		MappedRead * read = 0;
		if (!(readid & 0x1)) // Requested first read of a pair (even-numbered read id)
		{
			// Generate whole pair
			read = GenerateSingleRead(readid);
			nextRead = GenerateSingleRead(readid + 1);

			if (read != 0) {
				read->Paired = nextRead;
				if (nextRead == 0)
					read->SetFlag(NGMNames::NoSrcPair);
			}
			if (nextRead != 0) {
				nextRead->Paired = read;
				if (read == 0)
					nextRead->SetFlag(NGMNames::NoSrcPair);
			}
		} else {
			read = nextRead;
			nextRead = 0;

			if ((read != 0) && (read->ReadId != readid)) {
				Log.Error("Non-sequential read generation (req %i, got %i)", readid, read->ReadId);
				Fatal();
			}
		}
		return read;
	} else {
		return GenerateSingleRead(readid);
	}
}

void ReadProvider::DisposeRead(MappedRead * read) {
	if (NGM.Paired() && read->Paired != 0) // Program runs in paired mode and pair exists
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
