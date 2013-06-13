/*
 * main.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: philipp_
 */
#include "interleave-pairs.h"

#include <zlib.h>
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <map>
#include <cstring>
#include <tclap/CmdLine.h>

#include "Log.h"
#include "kseq.h"
#include "IParser.h"
#include "FastxParser.h"
#include "SamParser.h"
#include "Writer.h"
#include "FastqWriter.h"
#include "MappedRead.h"

#ifdef _BAM
#include "BamParser.h"
#endif

using std::cout;
using std::endl;
using std::map;
using std::string;

#undef module_name
#define module_name "MAIN"

struct nString {
	char * s;
	size_t n;
};

struct Sequence {
	nString name, seq, qlty;
	int fragPos;
};

char seperator = '/';

IParser * DetermineParser(char const * fileName) {
	IParser * parser = 0;
	gzFile fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("File %s does not exist!", fileName);
		Fatal();
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
		Log.Message("%s is SAM", fileName);
		parser = new SamParser();
	} else {
		if (strncmp(buffer, "BAM", 3) == 0) {
#ifdef _BAM
			Log.Message("%s is BAM", fileName);
			parser= new BamParser();
#else
			Log.Error("BAM input detected. NGM was compiled without BAM support!");
			Fatal();
#endif
		} else {
			if (buffer[0] == '>') {
				Log.Message("%s is FASTA", fileName);
			} else {
				Log.Message("%s is FASTQ", fileName);
			}
			parser = new FastXParser();
		}
	}
	delete[] buffer;
	return parser;
}

void cleareEntry(MappedRead * & entry) {
	if (entry != 0) {
//		if (entry->sequence != 0) {
//			delete entry->sequence;
//			entry->sequence = 0;
//		}
		delete entry;
		entry = 0;
	}
}

string parsePairId(string id) {
	return id.substr(0, id.rfind(seperator));
}

MappedRead * NextRead(IParser * parser, int const id) {

	MappedRead * read = 0;
	int l = parser->parseRead();
	//Log.Message("Name (%d): %s", seq->name.l, seq->name.s);
	//Log.Message("Seq  (%d): %s", seq->seq.l, seq->seq.s);
	//Log.Message("Qual (%d): %s", seq->qual.l, seq->qual.s);
	if (l >= 0) {
		if (parser->read->seq.l == parser->read->qual.l || parser->read->qual.l == 0) {
			read = new MappedRead(id, 10000);

			//Name
			static size_t const MAX_READNAME_LENGTH = 100;
			read->name = new char[MAX_READNAME_LENGTH];
			int nameLength = std::min(MAX_READNAME_LENGTH - 1, parser->read->name.l);
			memcpy(read->name, parser->read->name.s, nameLength);
			read->name[nameLength] = '\0';

//			char const * debugRead = "FCC01PDACXX:4:1101:10342:37018#0/1";
//			if(strcmp(read->name, debugRead) == 0) {
//				Log.Error("Read %s found: assigning id %d", debugRead, read->ReadId);
//			}

			//Sequence
			read->length = parser->read->seq.l;
			read->Seq = new char[read->length + 1];
			memset(read->Seq, '\0', read->length + 1);
			int nCount = 0;
			for (int i = 0; i < read->length; ++i) {
				char c = toupper(parser->read->seq.s[i]);
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
			if (parser->read->qual.l > 0) {
				read->qlty = new char[read->length + 1];
				memcpy(read->qlty, parser->read->qual.s, read->length);
				read->qlty[read->length] = '\0';
			}

//			Log.Message("%s", read->name);
//			Log.Message("%s", read->Seq);
//			if (read->qlty != 0)
//				Log.Message("%s", read->qlty);

		} else {
			Log.Error("Discarding read %s. Length of read not equal length of quality values.", parser->read->name.s);
			Fatal();
		}
	} else {
		if (l == -1) {
			Log.Verbose("End of input file reached.");
		} else {
			Log.Error("Error while parsing read %d (%d)", id, l);
			Fatal();
		}
	}
	return read;
}

int parseNext(map<string, MappedRead *> & pairs, MappedRead * read, Writer * writer, int fragNumber, int & writtenReads) {
//	Entry * entryLeft = new Entry;
//	entryLeft->first = isFirst;

	MappedRead * entryLeft = read;
	MappedRead * entryRight = 0;

	string nameLeft;
//
//	if (parser->read(nameLeft, entryLeft)) {
	if (read != 0) {
//		Entry * entryRight = 0;
		nameLeft = parsePairId(entryLeft->name);
		if (pairs.find(nameLeft) != pairs.end()) {
			Log.Verbose("Pair found %s", nameLeft.c_str());
			if (entryLeft->ReadId == 0) {
				entryRight = pairs[nameLeft];
				if (pairs[nameLeft]->ReadId == 0) {
					cout << "Warning: Found read pair in a single file (Name: " << nameLeft
					<< ")! Do you have two reads with the same name in one file?" << std::endl;
				}
			} else {
				entryRight = entryLeft;
				entryLeft = pairs[nameLeft];
				if (!pairs[nameLeft]->ReadId == 0) {
					cout << "Warning: Found read pair in a single file (Name: " << nameLeft
					<< ")! Do you have two reads with the same name in one file?" << std::endl;
				}
			}
//
			writer->writeRead(entryLeft);
			writer->writeRead(entryRight);
//			//writerLeft->writeRead(seq1->name.s, seq1->name.l, seq1->seq.s, seq1->seq.l, seq1->qual.s, seq1->qual.l);
//
//			//writerRight->writeRead(seq2->name.s, seq2->name.l, seq2->seq.s, seq2->seq.l, seq2->qual.s, seq2->qual.l);
//
//			//FilePosition outputPositionLeft = writerLeft->writeRead(entryLeft, maxReadLength);
//
//			//TODO: Ask if here should really be entryLeft???
//			//FilePosition outputPositionRight = writerRight->writeRead(entryRight, maxReadLength);
//
			writtenReads += 1;
////			if (outputPositionLeft != outputPositionRight) {
////				//TODO:
////				throw "Invalid output position in function parseNext.";
////			}
//			//writerLeft.writeIndexEntryRead(entryLeft, checkRead(entryLeft));
//			//writerRight.writeIndexEntryRead(entryRight, checkRead(entryRight));
//
			pairs.erase(nameLeft);
			cleareEntry(entryLeft);
			cleareEntry(entryRight);
			Log.Verbose("Pair deleted %s", nameLeft.c_str());
			return 1;
		} else {
			pairs[nameLeft] = read;
			Log.Verbose("Mate %d added: %s", read->ReadId, nameLeft.c_str());
		}
	} else {
		cleareEntry(entryLeft);
	}
	return 0;
}

int interleave_pairs(int argc, char **argv) {

	try {

		TCLAP::CmdLine cmd("Interleaves paired end reads from two FASTA/Q files into one FASTQ file.", ' ', "0.1", false);

		TCLAP::ValueArg<std::string> leftArg("1", "m1", "Upstream mates (FASTA/Q)", true, "", "file");
		TCLAP::ValueArg<std::string> rightArg("2", "m2", "Downstream mates (FASTA/Q)", true, "", "file");

		TCLAP::ValueArg<std::string> outArg("o", "output", "Output file", true, "", "file");
		TCLAP::ValueArg<std::string> unpairedArg("u", "unpaired", "Write reads without mate to this file.", false, "", "file");

		TCLAP::ValueArg<char> delimiterArg("d", "delimiter", "The character that precedes the 1 and 2 in the input files.", false, '/',
				"char");

		TCLAP::SwitchArg noprogressArg("", "noprogress", "Suppress progress output.", cmd, false);

		cmd.add(delimiterArg);
		cmd.add(unpairedArg);
		cmd.add(outArg);
		cmd.add(rightArg);
		cmd.add(leftArg);

		cmd.parse(argc, argv);

		//Log.Message("Interleave %s and %s to %s (delimiter %c)", leftArg.getValue().c_str(), rightArg.getValue().c_str(), outArg.getValue().c_str(), delimiterArg.getValue());


		_log = &Log;
		_Log::Init(); // Inits logging to file

		seperator = delimiterArg.getValue();

		Writer * writer = new FastqWriter(outArg.getValue().c_str());

		IParser * parser1 = DetermineParser(leftArg.getValue().c_str());
		parser1->init(leftArg.getValue().c_str());
		IParser * parser2 = DetermineParser(rightArg.getValue().c_str());
		parser2->init(rightArg.getValue().c_str());

		int l1 = 0;
		int l2 = 0;

		bool progess = !noprogressArg.getValue();

		map<string, MappedRead *> pairs;

		int pairNumber = 0;
		bool eof = false;
		MappedRead * read1 = 0;
		MappedRead * read2 = 0;
		int count = 0;
		int nPairs = 0;
		int nReads = 0;
		while (!eof) {
			read1 = NextRead(parser1, 0);
			if (read1 != 0)
				nReads++;
			read2 = NextRead(parser2, 1);
			if (read2 != 0)
				nReads++;

			nPairs += parseNext(pairs, read1, writer, 0, pairNumber);
			nPairs += parseNext(pairs, read2, writer, 1, pairNumber);

			count += 1;
			eof = (read1 == 0 && read2 == 0);

			if (count % 10000 == 0 && progess) {
				Log.Progress("Processed: %d", count);
				if(nPairs == 0) {
					Log.Error("0 proper pairs were found in the last 20000 reads. Please check if the seperator is set correctly or if you are using the correct files.");
					Fatal();
				}
			}
		}
		delete writer;
		writer = 0;
		int nUnmappedRead = 0;
		if (unpairedArg.getValue() != "" && pairs.size() > 0) {
			Log.Message("Writing unpaired reads to %s", unpairedArg.getValue().c_str());
			Writer * unmappedWriter = new FastqWriter(unpairedArg.getValue().c_str());
			for (map<string, MappedRead *>::iterator it = pairs.begin(); it != pairs.end(); it++) {
				MappedRead * read = it->second;
				unmappedWriter->writeRead(read);
				nUnmappedRead += 1;
			}
			delete unmappedWriter; unmappedWriter = 0;
		} else {
			nUnmappedRead = pairs.size();
		}

		Log.Message("Reads found in files: %d", nReads);
		Log.Message("Proper paires found: %d", nPairs);
		Log.Message("Unpaired reads found: %d", nUnmappedRead);

	} catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	} catch (std::ios_base::failure &e) {
		std::cerr << "Error: " << e.what() << std::endl;
	} catch(char const * msg) {
		std::cerr << "Exception: " << msg << std::endl;
	}

	return 0;
}
