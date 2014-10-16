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

int const maxReadLength = 1000;

IParser * DetermineParserStr(string strfileName) {
	char const * fileName = strfileName.c_str();
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
	gzclose(fp);

	int count = 0;
	for (size_t i = 0; i < 1000 && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
		if (buffer[i] == '\t') {
			count++;
		}
	}
	if (count >= 10) {
		Log.Message("%s is SAM", fileName);
		parser = new SamParser(maxReadLength);
	} else {
		if (strncmp(buffer, "BAM", 3) == 0) {
#ifdef _BAM
			Log.Message("%s is BAM", fileName);
			parser= new BamParser(maxReadLength);
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
			parser = new FastXParser(maxReadLength);
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
	MappedRead * read = new MappedRead(id, maxReadLength);

	int l = 0;

	try {
		l = parser->parseRead(read);

		if (l >= 0) {
			//Reduce memory usage by only
//			char * tmp = read->Seq;
//			read->Seq = new char[read->length + 1];
//			strcpy(read->Seq, tmp);
//			delete[] tmp; tmp = 0;

		} else {

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

int parseNext(map<string, MappedRead *> & pairs, MappedRead * read, Writer * writer, int fragNumber, int & writtenReads) {

	MappedRead * entryLeft = read;
	MappedRead * entryRight = 0;

	string nameLeft;
	if (read != 0) {

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

			writer->writeRead(entryLeft);
			writer->writeRead(entryRight);

			writtenReads += 1;

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

		TCLAP::SwitchArg forceArg("f", "force", "Force finishing even if no pairs are found.", cmd, false);

		cmd.add(delimiterArg);
		cmd.add(unpairedArg);
		cmd.add(outArg);
		cmd.add(rightArg);
		cmd.add(leftArg);

		cmd.parse(argc, argv);

		//Log.Message("Interleave %s and %s to %s (delimiter %c)", leftArg.getValue().c_str(), rightArg.getValue().c_str(), outArg.getValue().c_str(), delimiterArg.getValue());

		_log = &Log;
		_Log::Init(0, 0); // Inits logging to file

		seperator = delimiterArg.getValue();

		Writer * writer = new FastqWriter(outArg.getValue().c_str());

		IParser * parser1 = DetermineParserStr(leftArg.getValue());
		parser1->init(leftArg.getValue().c_str(), false);
		IParser * parser2 = DetermineParserStr(rightArg.getValue());
		parser2->init(rightArg.getValue().c_str(), false);

		int l1 = 0;
		int l2 = 0;

		bool progess = !noprogressArg.getValue();
		bool force = forceArg.getValue();

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
				if(nPairs == 0 && !force) {
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
				delete read; read = 0;
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
