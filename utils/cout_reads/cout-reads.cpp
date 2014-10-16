/*
 * cout-reads.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: fritz
 */




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
#include "SAMRecord.h"
#include "../paired/interleave-pairs.h"


#ifdef _BAM
#include "BamParser.h"
#endif

using std::cout;
using std::endl;
using std::map;
using std::string;

#undef module_name
#define module_name "MAIN"


//
//struct nString {
//	char * s;
//	size_t n;
//};
//
//struct Sequence {
//	nString name, seq, qlty;
//	int fragPos;
//};

int NextSAMRead(IParser * parser, int const id,SAMRecord * & read) {
//	cout<<"get next read"<<endl;

	int l = parser->parseSAMRecord(read);
	//Log.Message("Name (%d): %s", seq->name.l, seq->name.s);
	//Log.Message("Seq  (%d): %s", seq->seq.l, seq->seq.s);
	//Log.Message("Qual (%d): %s", seq->qual.l, seq->qual.s);
	if (l <= 0) {
		if (l == -2) {
			Log.Error("Read %s: Length of read not equal length of quality values.", read->get_read_name().c_str());
			Fatal();
		} else if (l != -1) {
			Log.Error("Unknown error while parsing read %d (%d)", id + 1, l);
			Fatal();
		}
		delete read;
		read = 0;
	}
//	cout<<"return"<<endl;
	return l;
}



int count_reads(int argc, char **argv) {
	try {

		TCLAP::CmdLine cmd(
				"Counts the number of mapped read per Chromosom using the sam/bam file",
				' ', "0.1", false);

		TCLAP::ValueArg<std::string> mapped_file("q", "s/bam",
				"Mapped read file (sam/bam", true, "", "file");

		TCLAP::ValueArg<std::string> outArg("o", "output", "Output file", true,
				"", "file");

//		TCLAP::SwitchArg noprogressArg("", "noprogress",
//				"Suppress progress output.", cmd, false);
		cmd.add(outArg);
		cmd.add(mapped_file);

//		cmd.add(noprogressArg);

		cmd.parse(argc, argv);

		_log = &Log;
		_Log::Init(0,0); // Inits logging to file


		IParser * parser1 = DetermineParserStr(mapped_file.getValue());

		parser1->init(mapped_file.getValue().c_str(), false);

		int l1 = 0;

		bool progess = true;//!noprogressArg.getValue();


		map<string, int> chrs;
		int nReads=0;
		bool eof = false;
		SAMRecord * read = new SAMRecord();

		while (NextSAMRead(parser1, 0,read) >= 0) {
			if (read != 0){
				nReads++;
				string chr=read->get_chr();
				if(chrs.find(chr)==chrs.end()){
					chrs[chr.c_str()]=1;
				}else{
					chrs[chr.c_str()]++;
				}
			}
			if (nReads % 10000 == 0 && progess) {
				Log.Progress("Processed: %d", nReads);
			}
		}
		FILE *file;
		file = fopen(outArg.getValue().c_str(), "w");
		//TODO: output file:
		for(map<string,int>::iterator i = chrs.begin();i!=chrs.end();i++){
			fprintf(file, "%s",(*i).first.c_str());
			fprintf(file, "%c",'\t');
			fprintf(file, "%i",(*i).second);
			fprintf(file, "%c",'\t');
			fprintf(file, "%f",(float)(*i).second/(float)nReads *100);
			fprintf(file, "%c",'\n');
		}
		fclose(file);
		Log.Message("Reads found in files: %d", nReads);

	} catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	} catch (std::ios_base::failure &e) {
		std::cerr << "Error: " << e.what() << std::endl;
	} catch(char const * msg) {
		std::cerr << "Exception: " << msg << std::endl;
	}

	return 0;
}
