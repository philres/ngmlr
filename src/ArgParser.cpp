#include "ArgParser.h"

#include <string.h>
#include <tclap/CmdLine.h>
#include <sstream>

#include "Log.h"
#include "PlatformSpecifics.h"

#undef module_name
#define module_name "ARGPARSER"

ArgParser::ArgParser(int argc, char * * argv) {
	ParseArguments(argc, (char const * *) argv);
}

ArgParser::~ArgParser() {

}

char * fromString(std::string str) {
	char * cstr = 0;
	if(str.size() > 0) {
		cstr = new char[str.size() + 1];
		strcpy(cstr, str.c_str());
	}
	return cstr;
}

void ArgParser::ParseArguments(int argc, char const * argv[]) {

	argv[0] = "ngmlr";

	TCLAP::CmdLine cmd("", ' ', "", false);

	TCLAP::ValueArg<std::string> queryArg("q", "query", "Path to the read file (FASTA/Q, SAM/BAM)", true, "", "file");
	TCLAP::ValueArg<std::string> refArg("r", "reference", "Path to the reference genome (FASTA/Q, can be gzipped)", true, "", "file");
	TCLAP::ValueArg<std::string> outArg("o", "output", "Path to output file (stdout if omitted)", false, "", "file");
	TCLAP::ValueArg<std::string> vcfArg("", "vcf", "SNPs will be taken into account when building reference index", false, "", "file");
	TCLAP::ValueArg<std::string> bedfilterArg("", "bed-filter", "If specified, only reads in the regions specified by the BED file are read from the input BAM file (requires BAM input)", false, "", "file");

	TCLAP::ValueArg<float> minIdentityArg("i", "min-identity", "All reads mapped with an identity lower than this threshold will be reported as unmapped (default: 0.65)", false, minIdentity, "0-1");
	TCLAP::ValueArg<float> minResiduesArg("R", "min-residues", "All reads mapped with lower than <int> or <float> * read length residues will be reported as unmapped. (default: 50)", false, minResidues, "int/float");
	TCLAP::ValueArg<float> sensitivityArg("s", "sensitivity", "", false, sensitivity, "0-1");

	TCLAP::ValueArg<int> threadsArg("t", "threads", "Number of threads", false, 1, "int");

	TCLAP::SwitchArg noprogressArg("", "no-progress", "Don't print progress info while mapping", cmd, false);

	cmd.add(queryArg);
	cmd.add(refArg);
	cmd.add(outArg);
	cmd.add(bedfilterArg);
	cmd.add(vcfArg);
	cmd.add(threadsArg);
	cmd.add(minIdentityArg);
	cmd.add(minResiduesArg);

	cmd.parse(argc, argv);

	queryFile = fromString(queryArg.getValue());
	referenceFile = fromString(refArg.getValue());
	outputFile = fromString(outArg.getValue());
	vcfFile = fromString(vcfArg.getValue());
	bedFile = fromString(bedfilterArg.getValue());

	minIdentity = minIdentityArg.getValue();
	minResidues = minResiduesArg.getValue();

	threads = threadsArg.getValue();
	progress = !noprogressArg.getValue();

	std::stringstream fullCmd;
	fullCmd << std::string(argv[0]);
	for(int i = 1; i < argc; ++i) {
		fullCmd << " " << std::string(argv[i]);
	}
	fullCommandLineCall = fromString(fullCmd.str());

	if(!FileExists(queryFile)) {
		Log.Error("Query file (%s) does not exist.", queryFile);
	}
	if(!FileExists(referenceFile)) {
		Log.Error("Reference file (%s) does not exist.", referenceFile);
	}
	if(bedFile != 0 && !FileExists(bedFile)) {
		Log.Error("BED filter file (%s) does not exist.", bedFile);
	}
	if(vcfFile != 0 && !FileExists(vcfFile)) {
		Log.Error("SNP file (%s) does not exist.", vcfFile);
	}
}

