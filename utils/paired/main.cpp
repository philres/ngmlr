/*
 * main.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: philipp_
 */
#include <zlib.h>
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <map>
#include <cstring>

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

ILog const * _log = 0;

struct nString {
	char * s;
	size_t n;
};

struct Sequence {
	nString name, seq, qlty;
	int fragPos;
};

char seperator = '/';

// actually platform specific.../care
ulong const FileSize(char const * const filename) {
	FILE * fp = fopen(filename, "rb");
	if (fp == 0) {
		Log.Warning("Tried to get size of nonexistant file %s", filename);
		return 0;
	}

	if (fseek(fp, 0, SEEK_END) != 0)
		return 0;

	ulong end = ftell(fp);
	fclose(fp);
	return end;
}

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

MappedRead * NextRead(IParser * parser, kseq_t *seq, int const id) {

	MappedRead * read = 0;
	int l = parser->parseRead(seq);
	//Log.Message("Name (%d): %s", seq->name.l, seq->name.s);
	//Log.Message("Seq  (%d): %s", seq->seq.l, seq->seq.s);
	//Log.Message("Qual (%d): %s", seq->qual.l, seq->qual.s);
	if (l >= 0) {
		if (seq->seq.l == seq->qual.l || seq->qual.l == 0) {
			read = new MappedRead(id, 10000);

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
			read->length = seq->seq.l;
			read->Seq = new char[read->length + 1];
			memset(read->Seq, '\0', read->length + 1);
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

int main(int argc, char **argv) {
	kseq_t *seq1;
	kseq_t *seq2;

	_log = &Log;
	_Log::Init(); // Inits logging to file

	Writer * writer = new FastqWriter(argv[3]);

	IParser * parser1 = DetermineParser(argv[1]);
	seq1 = parser1->init_seq(argv[1]);
	IParser * parser2 = DetermineParser(argv[2]);
	seq2 = parser2->init_seq(argv[2]);

	int l1 = 0;
	int l2 = 0;
//	int count = 0;
//	while ((l1 = parser1->parseRead(seq1)) > 0 && (l2 = parser2->parseRead(seq2)) > 0) {
//
//		//Log.Message("Name: %s", seq1->name.s);
//		//Log.Message("Read: %s", seq1->seq.s);
//		writer->writeRead(seq1->name.s, seq1->name.l, seq1->seq.s, seq1->seq.l, seq1->qual.s, seq1->qual.l);
//
//		//Log.Message("Name: %s", seq2->name.s);
//		//Log.Message("Read: %s", seq2->seq.s);
//		writer->writeRead(seq2->name.s, seq2->name.l, seq2->seq.s, seq2->seq.l, seq2->qual.s, seq2->qual.l);
//
//		if (count++ > 1000) {
//			return 0;
//		}
//	}

	map<string, MappedRead *> pairs;

	int pairNumber = 0;
	bool eof = false;
	MappedRead * read1 = 0;
	MappedRead * read2 = 0;
	int count = 0;
	int nPairs = 0;
	int nReads = 0;
	while (!eof) {
		read1 = NextRead(parser1, seq1, 0);
		if (read1 != 0)
			nReads++;
		read2 = NextRead(parser2, seq2, 1);
		if (read2 != 0)
			nReads++;

		nPairs += parseNext(pairs, read1, writer, 0, pairNumber);
		nPairs += parseNext(pairs, read2, writer, 1, pairNumber);

		count += 1;
		eof = (read1 == 0 && read2 == 0);

		if(count % 10000 == 0) {
			Log.Progress("Processed: %d", count);
		}
	}

	Log.Message("Reads read: %d", nReads);
	Log.Message("Proper paires found: %d", nPairs);

//		int writtenReads = 2 * pairNumber + pairs.size();
//		int lineNumber = pairNumber + pairs.size();
//		if (pairs.size() > 0) {
//			for (map<string, Entry *>::iterator it = pairs.begin(); it != pairs.end(); it++) {
//				Entry * entry = it->second;
//
//				bool valid = checkRead(entry);
//				if (entry->first) {
//					//writerLeft.writeIndexEntryRead(entry, valid);
//					writerLeft->writeRead(entry, maxReadLengthLeft);
//					//left
//					for (vector<Statistics*>::iterator i = jops_left.begin(); i != jops_left.end(); i++) {
//						(*(*i)).add(entry->sequence);
//					}
//					writerRight->writeDummy(maxReadLengthLeft);
//					//Dummy entry: has to be ignored
//					//writerRight.writeIndexEntryRead(entry, false);
//				} else {
//					writerRight->writeRead(entry, maxReadLengthRight);
//					//writerRight.writeIndexEntryRead(entry, valid);
//					//right
//					for (vector<Statistics*>::iterator i = jops_right.begin(); i != jops_right.end(); i++) {
//						(*(*i)).add(entry->sequence);
//					}
//					writerLeft->writeDummy(maxReadLengthRight);
//					//Dummy entry: has to be ignored
//					//writerLeft.writeIndexEntryRead(entry, false);
//				}
//				//			if (outputPosition != outputPositionDummy) {
//				//				//TODO:
//				//				throw "Invalid output position in function convertPairs.";
//				//			}
//			}
//		}

	return 0;
}
