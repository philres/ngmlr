/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "ArgParser.h"

#include <string.h>
#include <tclap/CmdLine.h>
#include <sstream>

#include "ArgParseOutput.h"
#include "Version.h"
#include "Log.h"
#include "PlatformSpecifics.h"

#undef module_name
#define module_name "ARGPARSER"

ArgParser::ArgParser(int argc, char * * argv) {

	outDefault = "stdout";
	noneDefault = "none";

	ParseArguments(argc, (char const * *) argv);
}

ArgParser::~ArgParser() {

}

char * ArgParser::fromString(std::string str) {
	char * cstr = 0;
	if(str.size() > 0 && str != outDefault && str != noneDefault) {
		cstr = new char[str.size() + 1];
		strcpy(cstr, str.c_str());
	}
	return cstr;
}

template<typename T>
void printParameter(std::stringstream & usage, TCLAP::ValueArg<T> & arg) {

	usage << "    " <<  arg.longID() << std::endl;
	usage << "        " << arg.getDescription();
	if(!arg.isRequired()) {
		usage << " [" << arg.getValue() << "]";
	}
	usage << std::endl;
}

void printParameter(std::stringstream & usage, TCLAP::SwitchArg & arg) {

	usage << "    " <<  arg.longID() << std::endl;
	usage << "        " << arg.getDescription();
	if(!arg.isRequired()) {
		usage << " [" << (arg.getValue() ? "true" : "false") << "]";
	}
	usage << std::endl;
}

void ArgParser::ParseArguments(int argc, char const * argv[]) {

	argv[0] = "ngmlr";

	TCLAP::CmdLine cmd("", ' ', "", true);

	TCLAP::ValueArg<std::string> queryArg("q", "query", "Path to the read file (FASTA/Q)", false, "/dev/stdin", "file", cmd);
	TCLAP::ValueArg<std::string> refArg("r", "reference", "Path to the reference genome (FASTA/Q, can be gzipped)", true, noneDefault, "file", cmd);
	TCLAP::ValueArg<std::string> outArg("o", "output", "Path to output file", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> vcfArg("", "vcf", "SNPs will be taken into account when building reference index", false, noneDefault, "file", cmd);
	TCLAP::ValueArg<std::string> bedfilterArg("", "bed-filter", "Only reads in the regions specified by the BED file are read from the input file (requires BAM input)", false, noneDefault, "file", cmd);

	TCLAP::ValueArg<std::string> rgIdArg("", "rg-id", "Adds RG:Z:<string> to all alignments in SAM/BAM", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgSmArg("", "rg-sm", "RG header: Sample", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgLbArg("", "rg-lb", "RG header: Library", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgPlArg("", "rg-pl", "RG header: Platform", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgDsArg("", "rg-ds", "RG header: Description", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgDtArg("", "rg-dt", "RG header: Date (format: YYYY-MM-DD)", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgPuArg("", "rg-pu", "RG header: Platform unit", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgPiArg("", "rg-pi", "RG header: Median insert size", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgPgArg("", "rg-pg", "RG header: Programs", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgCnArg("", "rg-cn", "RG header: sequencing center", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgFoArg("", "rg-fo", "RG header: Flow order", false, noneDefault, "string", cmd);
	TCLAP::ValueArg<std::string> rgKsArg("", "rg-ks", "RG header: Key sequence", false, noneDefault, "string", cmd);

	TCLAP::ValueArg<std::string> presetArgs("x", "presets", "Parameter presets for different sequencing technologies", false, "pacbio", "pacbio, ont", cmd);

	TCLAP::ValueArg<float> minIdentityArg("i", "min-identity", "Alignments with an identity lower than this threshold will be discarded", false, minIdentity, "0-1", cmd);
	TCLAP::ValueArg<float> minResiduesArg("R", "min-residues", "Alignments containing less than <int> or (<float> * read length) residues will be discarded", false, minResidues, "int/float", cmd);
	TCLAP::ValueArg<float> sensitivityArg("s", "sensitivity", "", false, sensitivity, "0-1", cmd);

	TCLAP::ValueArg<int> threadsArg("t", "threads", "Number of threads", false, 1, "int", cmd);

	TCLAP::ValueArg<int> binSizeArg("", "bin-size", "Sets the size of the grid used during candidate search", false, binSize, "int", cmd);
	TCLAP::ValueArg<int> kmerLengthArg("k", "kmer-length", "K-mer length in bases", false, kmerLength, "10-15", cmd);
	TCLAP::ValueArg<int> kmerSkipArg("", "kmer-skip", "Number of k-mers to skip when building the lookup table from the reference", false, kmerSkip, "int", cmd);
	TCLAP::ValueArg<int> maxInitialSegmentsArg("", "max-segments", "Max number of segments allowed for a read per kb", false, maxSegmentNumberPerKb, "int", cmd);
	TCLAP::ValueArg<float> scoreMatchArg("", "match", "Match score", false, scoreMatch, "float", cmd);
	TCLAP::ValueArg<float> scoreMismatchArg("", "mismatch", "Mismatch score", false, scoreMismatch, "float", cmd);
	TCLAP::ValueArg<float> scoreGapOpenArg("", "gap-open", "Gap open score", false, scoreGapOpen, "float", cmd);
	TCLAP::ValueArg<float> scoreGapExtendMaxArg("", "gap-extend-max", "Gap open extend max", false, scoreGapExtendMax, "float", cmd);
	TCLAP::ValueArg<float> scoreGapExtendMinArg("", "gap-extend-min", "Gap open extend min", false, scoreGapExtendMin, "float", cmd);
	TCLAP::ValueArg<float> scoreGapDecayArg("", "gap-decay", "Gap extend decay", false, scoreGapDecay, "float", cmd);
	TCLAP::ValueArg<int> stdoutArg("", "stdout", "Debug mode", false, stdoutMode, "0-7", cmd);
	TCLAP::ValueArg<int> subreadAligner("", "subread-aligner", "Choose subread aligning method", false, subreadaligner, "0-3", cmd);
	TCLAP::ValueArg<int> readpartLengthArg("", "subread-length", "Length of fragments reads are split into", false, readPartLength, "int", cmd);
	TCLAP::ValueArg<int> readpartCorridorArg("", "subread-corridor", "Length of corridor sub-reads are aligned with", false, readPartCorridor, "int", cmd);
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
	TCLAP::SwitchArg lowqualitysplitArg("", "no-lowqualitysplit", "Split alignments with poor quality", cmd, false);
	TCLAP::SwitchArg nosmallInversionArg("", "no-smallinv", "Don't detect small inversions", cmd, false);
	TCLAP::SwitchArg printAllAlignmentsArg("", "print-all", "Print all alignments. Disable filtering. (debug)", cmd, false);
	TCLAP::SwitchArg skipWriteArg("", "skip-write", "Don't write reference index to disk", cmd, skipSave);
	//skipSave = false;
	//updateCheck = false;
	//writeUnmapped = true;
	TCLAP::SwitchArg noSSEArg("", "nosse", "Debug switch (don't use if you don't know what you are doing)", cmd, false);
	TCLAP::SwitchArg bamFixArg("", "bam-fix", "Report reads with > 64k CIGAR operations as unmapped. Required to be compatibel to BAM format", cmd, false);
	TCLAP::SwitchArg skipAlignArg("", "skip-align", "Skip alignment step. Only for debugging purpose", cmd, false);

	std::stringstream usage;
	usage << "" << std::endl;
	usage << "Usage: ngmlr [options] -r <reference> -q <reads> [-o <output>]" << std::endl;
	usage << "" << std::endl;
	usage << "Input/Output:" << std::endl;
	printParameter<std::string>(usage, refArg);
	printParameter<std::string>(usage, queryArg);
	printParameter<std::string>(usage, outArg);
	printParameter(usage, skipWriteArg);
	printParameter(usage, bamFixArg);
	printParameter<std::string>(usage, rgIdArg);
	printParameter<std::string>(usage, rgSmArg);
	printParameter<std::string>(usage, rgLbArg);
	printParameter<std::string>(usage, rgPlArg);
	printParameter<std::string>(usage, rgDsArg);
	printParameter<std::string>(usage, rgDtArg);
	printParameter<std::string>(usage, rgPuArg);
	printParameter<std::string>(usage, rgPiArg);
	printParameter<std::string>(usage, rgPgArg);
	printParameter<std::string>(usage, rgCnArg);
	printParameter<std::string>(usage, rgFoArg);
	printParameter<std::string>(usage, rgKsArg);

	usage << "" << std::endl;
	usage << "General:" << std::endl;
	printParameter<int>(usage, threadsArg);
	printParameter<std::string>(usage, presetArgs);
	printParameter<float>(usage, minIdentityArg);
	printParameter<float>(usage, minResiduesArg);
	printParameter(usage, nosmallInversionArg);
	printParameter(usage, lowqualitysplitArg);
	printParameter(usage, verboseArg);
	printParameter(usage, noprogressArg);
	usage << "" << std::endl;
	usage << "Advanced:" << std::endl;
	printParameter<float>(usage, scoreMatchArg);
	printParameter<float>(usage, scoreMismatchArg);
	printParameter<float>(usage, scoreGapOpenArg);
	printParameter<float>(usage, scoreGapExtendMaxArg);
	printParameter<float>(usage, scoreGapExtendMinArg);
	printParameter<float>(usage, scoreGapDecayArg);
	printParameter<int>(usage, kmerLengthArg);
	printParameter<int>(usage, kmerSkipArg);
	printParameter<int>(usage, binSizeArg);
	printParameter<int>(usage, maxInitialSegmentsArg);
	printParameter<int>(usage, readpartLengthArg);
	printParameter<int>(usage, readpartCorridorArg);
	//printParameter<std::string>(usage, vcfArg);
	//printParameter<std::string>(usage, bedfilterArg);

	cmd.setOutput(new ArgParseOutput(usage.str(), ""));

	cmd.parse(argc, argv);

	queryFile = fromString(queryArg.getValue());
	referenceFile = fromString(refArg.getValue());
	outputFile = fromString(outArg.getValue());
	vcfFile = fromString(vcfArg.getValue());
	bedFile = fromString(bedfilterArg.getValue());
	rgId = fromString(rgIdArg.getValue());
	rgSm = fromString(rgSmArg.getValue());
	rgLb = fromString(rgLbArg.getValue());
	rgPl = fromString(rgPlArg.getValue());
	rgDs = fromString(rgDsArg.getValue());
	rgDt = fromString(rgDtArg.getValue());
	rgPu = fromString(rgPuArg.getValue());
	rgPi = fromString(rgPiArg.getValue());
	rgPg = fromString(rgPgArg.getValue());
	rgCn = fromString(rgCnArg.getValue());
	rgFo = fromString(rgFoArg.getValue());
	rgKs = fromString(rgKsArg.getValue());

	minIdentity = minIdentityArg.getValue();
	minResidues = minResiduesArg.getValue();

	binSize = binSizeArg.getValue();
	kmerLength = kmerLengthArg.getValue();
	threads = threadsArg.getValue();
	kmerSkip = kmerSkipArg.getValue();
	maxSegmentNumberPerKb = maxInitialSegmentsArg.getValue();
	scoreMatch = scoreMatchArg.getValue();
	if(scoreMatch < 0.0f) {
		Log.Message("--match must not be smaller than zero. changing from %f to %f", scoreMatch, scoreMatch * -1.0f);
		scoreMatch = scoreMatch * -1.0f;
	}
	scoreMismatch = scoreMismatchArg.getValue();
	if (scoreMismatch > 0.0f) {
		Log.Message("--mismatch must not be greater than zero. changing from %f to %f", scoreMismatch, scoreMismatch * -1.0f);
		scoreMismatch = scoreMismatch * -1.0f;
	}
	scoreGapOpen = scoreGapOpenArg.getValue();
	if (scoreGapOpen > 0.0f) {
		Log.Message("--gap-open must not be greater than zero. changing from %f to %f", scoreGapOpen, scoreGapOpen * -1.0f);
		scoreGapOpen = scoreGapOpen * -1.0f;
	}
	scoreGapExtendMax = scoreGapExtendMaxArg.getValue();
	if (scoreGapExtendMax > 0.0f) {
		Log.Message("--gap-extend-max must not be greater than zero. changing from %f to %f", scoreGapExtendMax, scoreGapExtendMax * -1.0f);
		scoreGapExtendMax = scoreGapExtendMax * -1.0f;
	}
	scoreGapExtendMin = scoreGapExtendMinArg.getValue();
	if (scoreGapExtendMin > 0.0f) {
		Log.Message("--gap-extend-min must not be greater than zero. changing from %f to %f", scoreGapExtendMin, scoreGapExtendMin * -1.0f);
		scoreGapExtendMin = scoreGapExtendMin * -1.0f;
	}
	scoreGapDecay = scoreGapDecayArg.getValue();
	if (scoreGapDecay < 0.0f) {
		Log.Message("--gap-decay must not be smaller than zero. changing from %f to %f", scoreGapDecay, scoreGapDecay * -1.0f);
		scoreGapDecay = scoreGapDecay * -1.0f;
	}
	subreadaligner = subreadAligner.getValue();
	stdoutMode = stdoutArg.getValue();
	readPartCorridor = readpartCorridorArg.getValue();
	readPartLength = readpartLengthArg.getValue();

	progress = !noprogressArg.getValue();
	color = colorArg.getValue();
	verbose = verboseArg.getValue();
	lowQualitySplit = !lowqualitysplitArg.getValue();
	smallInversionDetection = !nosmallInversionArg.getValue();
	printAllAlignments = printAllAlignmentsArg.getValue();
	skipSave = skipWriteArg.getValue();
	nosse = noSSEArg.getValue();
	bamCigarFix = bamFixArg.getValue();
	skipAlign = skipAlignArg.getValue();

	if (presetArgs.getValue() == "pacbio") {
		//Do nothing. Defaults are for Pacbio
	} else if (presetArgs.getValue() == "ont") {
//		lowQualitySplit = (lowqualitysplitArg.isSet()) ? lowQualitySplit : false;
//		smallInversionDetection = (nosmallInversionArg.isSet()) ? smallInversionDetection : false;
//		scoreMatch = (scoreMatchArg.isSet()) ? scoreMatch : 3;
//		scoreMismatch = (scoreMatchArg.isSet()) ? scoreMismatch : -3;
		//scoreGapOpen = (scoreGapOpenArg.isSet()) ? scoreGapOpen : -1;
		//scoreGapExtendMax = (scoreGapExtendMaxArg.isSet()) ? scoreGapExtendMax : -1;
		//scoreGapExtendMax = (scoreGapExtendMinArg.isSet()) ? scoreGapExtendMin : -0.5;
		scoreGapDecay = (scoreGapDecayArg.isSet()) ? scoreGapDecay : 0.15;
	} else {
		std::cerr << "Preset " << presetArgs.getValue() << " not found" << std::endl;
	}

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

