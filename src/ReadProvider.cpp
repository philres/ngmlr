/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "ReadProvider.h"

#include <zlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <limits.h>

#include "IConfig.h"
#include "Log.h"
#include "Timing.h"
#include "CS.h"
#include "MappedRead.h"
#include "FastxParser.h"
//#include "BamParser.h"
#include "SamParser.h"

using NGMNames::ReadStatus;

#undef module_name
#define module_name "INPUT"

ReadProvider::ReadProvider() :
		//bufferLength based on max read length of 1MB and read part length of readPartLength
		readPartLength(Config.getReadPartLength()), bufferLength(2000), parsedReads(0), /*readBuffer(
		 new MappedRead *[bufferLength]),*/readsInBuffer(0), parser1(0) {

}

uint ReadProvider::init() {
	typedef void (*PrefixIterationFn)(ulong prefix, uloc pos, ulong mutateFrom,
			ulong mutateTo, void* data);

	char const * const fileName1 = Config.getQueryFile();

	Log.Verbose("Initializing ReadProvider");

	Log.Message("Opening query file %s", fileName1);

	size_t maxLen = readPartLength + 16;
	parser1 = DetermineParser(fileName1, maxLen);

	return 0;
}

ReadProvider::~ReadProvider() {
	if (parser1 != 0) {
		delete parser1;
		parser1 = 0;
	}
}

void ReadProvider::splitRead(MappedRead * read) {


	int splitNumber = read->length / readPartLength;

//	int splitNumber = read->length / (readPartLength / 2);

	int nameLength = strlen(read->name);

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

		if ((readPartLength + 16) <= length) {
			Log.Message("ERror ereroeroeroeor");
		}
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

	} else {
		group->readNumber = splitNumber;
		group->reads = new MappedRead *[splitNumber];
		memset(group->reads, 0, sizeof(MappedRead *) * splitNumber);

		for (int i = splitNumber - 1; i >= 0; --i) {
			MappedRead * readPart = new MappedRead(read->ReadId + i,
					readPartLength + 16);

			strcpy(readPart->name, read->name);

			int length = std::min(readPartLength,
					read->length - i * readPartLength);
			readPart->length = length;

			readPart->Seq = new char[readPartLength + 16];
			memset(readPart->Seq, '\0', readPartLength + 16);
			strncpy(readPart->Seq, read->Seq + i * readPartLength, length);

//			strncpy(readPart->Seq, read->Seq + i * (readPartLength / 2), length);

			readPart->qlty = 0;

			readPart->group = group;

			readPart->group->reads[i] = readPart;

		}
	}
}

MappedRead * ReadProvider::NextRead(IParser * parser, int const id) {

	int l = 0;

//	if (readsInBuffer == 0 && ++parsedReads <= 20) {
//	if (readsInBuffer == 0) {
	try {
		static int const qryMaxLen = readPartLength + 16;
		MappedRead * read = new MappedRead(id, qryMaxLen);
		l = parser->parseRead(read);
//		Log.Message("Parsing next read: %s (%d)", read->name, read->length);
		if (l >= 0) {

			Log.Debug(2, "READ_%d\tINPUT\t%s", id, read->name);
			Log.Debug(16384, "READ_%d\tINPUT_DETAILS\t%s\t%s\t%s\t%s", id, read->Seq, read->qlty, read->AdditionalInfo);

			if(l > readPartLength) {
				splitRead(read);
			} else {
				read->group = 0;
			}

			NGM.AddReadRead(read->ReadId);

			NGM.Stats->readsInProcess += 1;
			return read;

		} else {

			Log.Debug(2, "READ_%d\tINPUT\t%s error while reading", id, read->name);

			if(l == -2) {
				Log.Error("Read %s: Length of read not equal length of quality values.", read->name);
			} else if (l != -1) {
				//TODO correct number when paired
				Log.Error("Unknown error while parsing read number %d (error code: %d)", id + 1, l);
			}
		}
		delete read;
		read = 0;
	} catch (char * ex) {
		Log.Error("%s", ex);
	}
//	}

//	if (readsInBuffer == 0) {
	return 0;
//	} else {
////		Log.Message("Sending already paresed read %d", readBuffer[readsInBuffer - 1]->ReadId);
//		return readBuffer[readsInBuffer-- - 1];
//	}
}

IParser * ReadProvider::DetermineParser(char const * fileName,
		int const qryMaxLen) {

	IParser * parser = 0;

	parser = new FastXParser(qryMaxLen);
	parser->init(fileName);

//	gzFile fp = gzopen(fileName, "r");
//	if (fp == 0) {
//		//File does not exist
//		Log.Error("File %s does not exist!",fileName);
//		throw "File not found.";
//	} else {
//		char * buffer = new char[1000];
//		while (gzgets(fp, buffer, 1000) > 0 && buffer[0] == '@') {
//		}
//
//		int count = 0;
//		for (size_t i = 0; i < 1000 && buffer[i] != '\0' && buffer[i] != '\n';
//				i++) {
//			if (buffer[i] == '\t') {
//				count++;
//			}
//		}
//		if (count >= 10) {
//			Log.Message("Input is SAM");
//			parser = new SamParser(qryMaxLen);
//		} else {
//			if (strncmp(buffer, "BAM", 3) == 0) {
////			Log.Message("Input is BAM");
////			parser= new BamParser(qryMaxLen);
//				Log.Error("BAM input is currently not supported!");
//			} else {
//				if (buffer[0] == '>') {
//					Log.Message("Input is Fasta");
//				} else {
//					Log.Message("Input is Fastq");
//				}
//				parser = new FastXParser(qryMaxLen);
//			}
//		}
//		gzclose(fp);
//		delete[] buffer;
//		buffer = 0;
//		parser->init(fileName);
//	}
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

//	Log.Message("Readseq: %s", read1->Seq);
	return read1 != 0;
}

void ReadProvider::DisposeRead(MappedRead * read) {
//	Log.Message("Disposing read %s", read->name);
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

	NGM.Stats->readsInProcess -= 1;

}
