/*
 * filter.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: philipp_
 */

#include "filter.h"

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
#include "../cout_reads/cout-reads.h"

#ifdef _BAM
#include "BamParser.h"
#endif

using std::cout;
using std::endl;
using std::map;
using std::string;

#undef module_name
#define module_name "MAIN"

int filter(int argc, char **argv) {
	try {

		TCLAP::CmdLine cmd(
				"Filter sam/bam file for % identity in the alignments",
				' ', "0.1", false);

		TCLAP::ValueArg<std::string> mapped_file("q", "s/bam",
				"Mapped read file (sam/bam", true, "", "file");

		TCLAP::ValueArg<float> identArg("i", "ident",
					"Min identity (0-100) %", true, "", "Number");

		TCLAP::ValueArg<int> MQArg("m", "mq",
							"Min mapping quality (0-256) %", true, "", "Number");

		TCLAP::ValueArg<std::string> outArg("o", "output", "Output sam/bam file", true,
				"", "file");

		cmd.add(MQArg);
		cmd.add(identArg);
		cmd.add(outArg);
		cmd.add(mapped_file);

		cmd.parse(argc, argv);

		_log = &Log;
		_Log::Init(); // Inits logging to file

		IParser * parser1 = DetermineParser(mapped_file.getValue().c_str());

		parser1->init(mapped_file.getValue().c_str(), false);

		int l1 = 0;

		bool progess = true; //!noprogressArg.getValue();

		map<string, int> chrs;
		int nReads = 0;
		bool eof = false;
		SAMRecord * read = new SAMRecord();

		while (NextSAMRead(parser1, 0, read) >= 0) {
			if (read != 0) {
				nReads++;
			}
			string chr = read->get_chr();
			if (chrs.find(chr) == chrs.end()) {
				chrs[chr.c_str()] = 1;
			} else {
				chrs[chr.c_str()]++;
			}

			if (nReads % 10000 == 0 && progess) {
				Log.Progress("Processed: %d", nReads);
			}
		}
		FILE *file;
		file = fopen(outArg.getValue().c_str(), "w");
		//TODO: output file:
		for (map<string, int>::iterator i = chrs.begin(); i != chrs.end(); i++) {
			fprintf(file, "%s", (*i).first.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", (*i).second);
			fprintf(file, "%c", '\t');
			fprintf(file, "%f", (float) (*i).second / (float) nReads * 100);
			fprintf(file, "%c", '\n');
		}
		fclose(file);
		Log.Message("Reads found in files: %d", nReads);

	}
	catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	} catch (std::ios_base::failure &e) {
		std::cerr << "Error: " << e.what() << std::endl;
	} catch(char const * msg) {
		std::cerr << "Exception: " << msg << std::endl;
	}

	return 0;
}

