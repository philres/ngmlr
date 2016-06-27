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
		//bufferLength based on max read length of 1MB and read part length of readPartLength
		readPartLength(Config.GetInt(READ_PART_LENGTH)), bufferLength(2000), parsedReads(0), readBuffer(
				new MappedRead *[bufferLength]), readsInBuffer(0), parser1(0), parser2(
				0), peDelimiter(
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
	size_t maxLen = readPartLength;
	bool estimate = !(Config.Exists("skip_estimate")
			&& Config.GetInt("skip_estimate"));
	if (!Config.Exists("qry_max_len") || estimate) {

		int avgLen = readPartLength;

//		if (maxLen > maxReadLength) {
//			Log.Warning("Max. supported read length is 1000bp. All reads longer than 1000bp will be hard clipped in the output.");
//			maxLen = maxReadLength;
//		}

		((_Config*) _config)->Override("qry_max_len", (int) maxLen);
		((_Config*) _config)->Override("qry_avg_len", (int) avgLen);

		Log.Message("Read part length: %d", readPartLength);

		if (Config.Exists("corridor")) {
			Log.Verbose("Corridor witdh overwritten!");
		} else {
			((_Config*) _config)->Default("corridor", (int) (5 + avgLen * 0.15));
		}
		Log.Message("Corridor width: %d", Config.GetInt("corridor"));

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
	//Log.Message("Estimating parameter took %.3fs", tmr.ET());

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

	if (readBuffer != 0) {
		delete[] readBuffer;
		readBuffer = 0;
	}
}

void ReadProvider::splitRead(MappedRead * read, int const qryMaxLen) {

	read->qryMaxLen = read->length + 1;
//	Log.Message("New read with length: %d", read->length);

	int splitNumber = read->length / readPartLength;

	int nameLength = strlen(read->name);

//	printf("Name: %s\n", read->name);
//	printf("ID: %d\n", read->ReadId);
//	printf("Length: %d\n", read->length);
//	printf("Parts: %d\n", splitNumber);
//	printf("Seq: %s\n", read->Seq);

	ReadGroup * group = new ReadGroup();
	group->fullRead = read;
	group->readId = read->ReadId;
	group->bestScoreSum = 0;
	group->fwdMapped = 0;
	group->reverseMapped = 0;
	group->readsFinished = 0;

	read->group = group;

	if (splitNumber == 0) {
		splitNumber = 1;
		group->readNumber = splitNumber;
		group->reads = new MappedRead *[splitNumber];
		MappedRead * readPart = new MappedRead(read->ReadId + 1,
				readPartLength + 16);

//		readPart->name = new char[nameLength + 1];
		strcpy(readPart->name, read->name);

		int length = read->length;
		readPart->length = length;

		readPart->Seq = new char[readPartLength + 16];
		memset(readPart->Seq, '\0', readPartLength + 16);
//		memset(readPart->Seq, 'N', readPartLength);
		strncpy(readPart->Seq, read->Seq, length);

		//readPart->qlty = new char[readPartLength + 1];
		//memset(readPart->qlty, '\0', readPartLength + 1);
		//strncpy(readPart->qlty, read->qlty + i * readPartLength, length);
		readPart->qlty = 0;

		readPart->group = group;
		group->reads[0] = readPart;

		readBuffer[readsInBuffer++] = readPart;
	} else {
		group->readNumber = splitNumber;
		group->reads = new MappedRead *[splitNumber];
		memset(group->reads, 0, sizeof(MappedRead *) * splitNumber);

		for (int i = splitNumber - 1; i >= 0; --i) {
			MappedRead * readPart = new MappedRead(read->ReadId + i,
					readPartLength + 16);

//			readPart->name = new char[nameLength + 1];
			strcpy(readPart->name, read->name);

			int length = std::min(readPartLength,
					read->length - i * readPartLength);
			readPart->length = length;

			readPart->Seq = new char[readPartLength + 16];
			memset(readPart->Seq, '\0', readPartLength + 16);
//			memset(readPart->Seq, 'N', readPartLength);
			strncpy(readPart->Seq, read->Seq + i * readPartLength, length);

			//readPart->qlty = new char[readPartLength + 1];
			//memset(readPart->qlty, '\0', readPartLength + 1);
			//strncpy(readPart->qlty, read->qlty + i * readPartLength, length);
			readPart->qlty = 0;

			readPart->group = group;

			readPart->group->reads[i] = readPart;

//		printf("Name: %s\n", readPart->name);
//		printf("ID: %d\n", readPart->ReadId);
//		printf("Length: %d\n", readPart->length);
//		printf("Seq: %s\n", readPart->Seq);

			readBuffer[readsInBuffer++] = readPart;
		}
	}
}

MappedRead * ReadProvider::NextRead(IParser * parser, int const id) {

	int l = 0;

//	if (readsInBuffer == 0 && ++parsedReads <= 200) {
	if (readsInBuffer == 0) {
		try {
			static int const qryMaxLen = Config.GetInt("qry_max_len");
			MappedRead * read = new MappedRead(id, qryMaxLen);
			l = parser->parseRead(read);
//			Log.Message("Parsing next read: %s (%d)", read->name, read->ReadId);
			if (l >= 0) {

				Log.Debug(2, "READ_%d\tINPUT\t%s", id, read->name);
				Log.Debug(16384, "READ_%d\tINPUT_DETAILS\t%s\t%s\t%s\t%s", id, read->Seq, read->qlty, read->AdditionalInfo);

				splitRead(read, qryMaxLen);
				NGM.AddReadRead(read->ReadId);

				NGM.Stats->readsInProcess += 1;

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
			}
			//delete read;
			//read = 0;
		} catch (char * ex) {
			Log.Error("%s", ex);
			Fatal();
		}
	}

	if (readsInBuffer == 0) {
		return 0;
	} else {
//		Log.Message("Sending already paresed read %d", readBuffer[readsInBuffer - 1]->ReadId);
		return readBuffer[readsInBuffer-- - 1];
	}
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

	read1 = GenerateSingleRead(readid1);

	return read1 != 0;
}

void ReadProvider::DisposeRead(MappedRead * read) {
	Timer tmr;
	tmr.ST();
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
		Log.Verbose("Disposing read %s", read->name);
		if (read->group != 0) {
			for (int j = 0; j < read->group->readNumber; ++j) {
				if (read->group->reads[j] != 0) {
					delete read->group->reads[j];
					read->group->reads[j] = 0;
				}
			}
			delete[] read->group->reads;
			read->group->reads = 0;
			Log.Verbose("Deleting group");
			delete read->group;
			read->group = 0;
		}

		// Single mode or no existing pair
		delete read;
	}
	NGM.Stats->readsInProcess -= 1;
	Log.Verbose("Deleting read took %fs", tmr.ET());
}
