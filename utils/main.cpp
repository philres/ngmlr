/*
 * main.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: philipp_
 */
#include <iostream>
#include <cstring>
#include <vector>
#include <tclap/CmdLine.h>

#include "Log.h"
#include "paired/interleave-pairs.h"
#include "kmers/kmer-distribution.h"

using std::string;
using TCLAP::ValuesConstraint;
using TCLAP::UnlabeledValueArg;

#undef module_name
#define module_name "MAIN"

#ifdef NDEBUG
bool cDebug = false;
#else
bool cDebug = true;
#endif

ILog const * _log = 0;

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

int main(int argc, char **argv) {

	string name = "ngm-utils";

	char const * arch[] = { "x86", "x64" };
	char const * build = (cDebug) ? " (DEBUG)" : "";
	Log.Message("ngm-utils: %s%s (build %s %s)", arch[sizeof(void*) / 4 - 1], build, __DATE__, __TIME__);

	try {

		TCLAP::CmdLine cmd("", ' ', "0.1", false);

		std::vector<std::string> allowed;
		allowed.push_back("interleave");
		allowed.push_back("filter");
		allowed.push_back("kmer");
		ValuesConstraint<string> allowedVals(allowed);

		UnlabeledValueArg<string> nolabel("program",
				"Name of the program you want to use. Available programs: interleave",
				true, "string", "name", &allowedVals);

		cmd.add(nolabel);

		//Change program name
		argv[0] = const_cast<char *>(name.c_str());
		cmd.parse(std::min(argc, 2), argv);

		if (nolabel.getValue() == allowed[0]) {
			char const * name = "ngm-utils interleave";
			argv[1] = const_cast<char *>(name);
			interleave_pairs(argc - 1, argv + 1);
		} else if (nolabel.getValue() == allowed[1]) {
			char const * name = "ngm-utils filter";
			argv[1] = const_cast<char *>(name);
			//filter(argc - 1, argv + 1);
		} else if (nolabel.getValue() == allowed[2]) {
			char const * name = "ngm-utils kmer";
			argv[1] = const_cast<char *>(name);
			kmer_distribution(argc - 1, argv + 1);
		} else {
			throw TCLAP::ArgException("Invalid value found", "program");
		}

	} catch (TCLAP::ArgException &e) // catch any exceptions
	{
		std::cerr << "Error: " << e.error() << " for arg " << e.argId()
				<< std::endl;
	}

	return 0;
}
