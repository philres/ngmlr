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
	for (size_t i = 0; i < 1000 && buffer[i] != '\0' && buffer[i] != '\n';
			i++) {
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


SAMRecord * NextRead(IParser * parser, int const id) {

	SAMRecord * read = 0;
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

	return read;
}



int count_reads(int argc, char **argv) {
	try {
		TCLAP::CmdLine cmd(
				"Interleaves paired end reads from two FASTA/Q files into one FASTQ file.",
				' ', "0.1", false);

		TCLAP::ValueArg<std::string> mapped_file("q", "s/bam",
				"Mapped read file (sam/bam", true, "", "file");

		TCLAP::ValueArg<std::string> outArg("o", "output", "Output file", true,
				"", "file");

		TCLAP::SwitchArg noprogressArg("", "noprogress",
				"Suppress progress output.", cmd, false);

		cmd.add(mapped_file);
		cmd.add(outArg);
		cmd.add(noprogressArg);

		cmd.parse(argc, argv);

		_log = &Log;
		_Log::Init(); // Inits logging to file


		IParser * parser1 = DetermineParser(mapped_file.getValue().c_str());
		parser1->init(mapped_file.getValue().c_str(), false);

		int l1 = 0;

		bool progess = !noprogressArg.getValue();


		map<string, int> chrs;
		int nReads=0;
		bool eof = false;
		SAMRecord * read1 = 0;

		int count = 0;

		while (!eof) {
			read1 = NextRead(parser1, 0); //-1??
			if (read1 != 0){
				nReads++;
			}
			string chr=read1->get_chr();
			if(chrs.find(chr)==chrs.end()){
				chrs[chr.c_str()]=1;
			}else{
				chrs[chr.c_str()]++;
			}
			count += 1;
			eof = (read1 == 0);

			if (count % 10000 == 0 && progess) {
				Log.Progress("Processed: %d", count);
			}
		}
		FILE *file;
		file = fopen(outArg.getValue().c_str(), "w");
		//TODO: output file:
		for(map<string,int>::iterator i = chrs.begin();i!=chrs.end();i++){
			fprintf(file, "%s",(*i).first.c_str());
			fprintf(file, "%c",'\t');
			fprintf(file, "%i",(*i).second);
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
