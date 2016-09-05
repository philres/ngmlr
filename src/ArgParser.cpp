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

	TCLAP::ValueArg<int> binSizeArg("", "bin-size", "Sets the size of the grid NextgenMap uses during CMR search to: 2^n (default: 4)", false, binSize, "int");
	TCLAP::ValueArg<int> kmerLengthArg("k", "kmer-length", "Kmer length in bases. (default: 13)", false, kmerLength, "10-15");
	TCLAP::ValueArg<int> kmerSkipArg("", "kmer-skip", "Number of k-mers to skip when building the lookup table from the reference (default: 2)", false, kmerSkip, "int");
	TCLAP::ValueArg<int> scoreMatchArg("", "match", "Match score", false, scoreMatch, "int");
	TCLAP::ValueArg<int> scoreMismatchArg("", "mismatch", "Mismatch score", false, scoreMismatch, "int");
	TCLAP::ValueArg<int> scoreGapOpenArg("", "gap-open", "Gap open score", false, scoreGapOpen, "int");
	TCLAP::ValueArg<int> scoreGapExtendArg("", "gap-extend", "Gap open extend", false, scoreGapExtend, "int");
	TCLAP::ValueArg<int> stdoutArg("", "stdout", "Debug mode", false, stdoutMode, "0-7");
	TCLAP::ValueArg<int> readpartLengthArg("", "subread-length", "Length of fragments reads are split to", false, readPartLength, "int");
	TCLAP::ValueArg<int> readpartCorridorArg("", "subread-corridor", "Length of corridor sub-reads are aligned with", false, readPartCorridor, "int");
	//csSearchTableLength = 0;
	//logLevel = 0; //16383, 255
	//minKmerHits = 0;
	//maxCMRs = INT_MAX;

	TCLAP::SwitchArg noprogressArg("", "no-progress", "Don't print progress info while mapping", cmd, false);
	TCLAP::SwitchArg verboseArg("", "verbose", "Debug output", cmd, false);
	//bam = false;
	TCLAP::SwitchArg colorArg("", "color", "Colored command line output", cmd, false);
	//hardClip = false;
	//log = false;
	TCLAP::SwitchArg nolowqualitysplitArg("", "no-lowqualitysplit", "Don't split alignments with poor quality", cmd, false);
	TCLAP::SwitchArg nosmallInversionArg("", "no-smallinv", "Don't detect small inversions", cmd, false);
	TCLAP::SwitchArg printAllAlignmentsArg("", "print-all", "Print all alignments. Disable filtering. (debug)", cmd, false);
	//skipSave = false;
	//updateCheck = false;
	//writeUnmapped = true;
	TCLAP::SwitchArg fastArg("", "fast", "Debug switch (don't use if you don't know what you do)", cmd, false);

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

	binSize = binSizeArg.getValue();
	kmerLength = kmerLengthArg.getValue();
	threads = threadsArg.getValue();
	kmerSkip = kmerSkipArg.getValue();
	scoreMatch = scoreMatchArg.getValue();
	scoreMismatch = scoreMismatchArg.getValue();
	scoreGapOpen = scoreGapOpenArg.getValue();
	scoreGapExtend = scoreGapExtendArg.getValue();
	stdoutMode = stdoutArg.getValue();
	readPartCorridor = readpartCorridorArg.getValue();
	readPartLength = readpartLengthArg.getValue();

	progress = !noprogressArg.getValue();
	color = colorArg.getValue();
	verbose = verboseArg.getValue();
	lowQualitySplit = !nolowqualitysplitArg.getValue();
	smallInversionDetection = !nosmallInversionArg.getValue();
	printAllAlignments = printAllAlignmentsArg.getValue();
	fast = fastArg.getValue();

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

