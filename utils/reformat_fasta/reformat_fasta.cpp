/*
 * reformat_fasta.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: fritz
 */

#include <zlib.h>
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <fstream>

#include <tclap/CmdLine.h>
#include "Log.h"

#include "../paired/interleave-pairs.h"

using namespace std;

int next_sequence(ifstream & myfile, string & name, string & seq) {
	//zlib
	char buf=' ';
	if(name.empty()){
		buf = myfile.get();
	}
	name.clear();
	seq.clear();
//	cout << "first: " << buf << endl;
	//parse header:
	while (myfile.good() && (buf != '\0' && buf != '\n')) {

		if (buf != '\0' && buf != ' '  && buf != '>') {
			name += buf;
		}
		buf = myfile.get();
	}
//parse seq:
	while (myfile.good())     // loop while extraction from file is possible
	{
		if (buf == '>') {     //header!
			return seq.size();
		} else if (buf != '\0' && buf != '\n' && buf != ' ') {
			seq += buf;
		}
		buf = myfile.get();
	}

	return seq.size();
}

int reformat_fasta(int argc, char **argv) {

	try {

		TCLAP::CmdLine cmd(
				"Reformats a fasta file to fit the requirements for NGM or other tools",
				' ', "0.1", false);

		TCLAP::ValueArg<std::string> fasta_file("f", "fasta",
				"Fasta input file ", true, "", "file");

		TCLAP::ValueArg<std::string> outArg("o", "output", "Output file", true,
				"", "file");

		cmd.add(outArg);
		cmd.add(fasta_file);

//		cmd.add(noprogressArg);

		cmd.parse(argc, argv);

		_log = &Log;
		_Log::Init(0, 0); // Inits logging to file

		ifstream myfile; //test2
		myfile.open(fasta_file.getValue().c_str(), ifstream::in);
		if (!myfile.good()) {
			cout << "No such file " << fasta_file.getValue().c_str() << endl;
			exit(0);
		}

		FILE *file;
		file = fopen(outArg.getValue().c_str(), "w");

		int l1 = 0;

		bool progess = true; //!noprogressArg.getValue();

		string name;
		string seq;
		int nReads = 0;

		while (next_sequence(myfile, name, seq) > 0) {
			if (!seq.empty()) {
				nReads++;
				fprintf(file, "%c", '>');
				fprintf(file, "%s", name.c_str());
				fprintf(file, "%c", '\n');
				for (size_t i = 1; i < seq.size() + 1; i++) {
					fprintf(file, "%c", seq[i - 1]);
					if (i % 200 == 0) {
						fprintf(file, "%c", '\n');
					}
				}
				if (seq.size() + 1 % 200 != 0) {
					fprintf(file, "%c", '\n');
				}
			}
			if (nReads % 10000 == 0 && progess) {
				Log.Progress("Processed: %d", nReads);
			}
		}
		myfile.close();
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
