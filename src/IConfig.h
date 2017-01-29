/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __ICONF_H__
#define __ICONF_H__

#include <climits>

class IConfig {

protected:

	/**
	 ******************
	 * Default values *
	 ******************
	 */

	float minIdentity = 0.65f;
	float minResidues = 0.25f;
	float sensitivity = 0.8f;

	int binSize = 4;
	int csSearchTableLength = 0;
	int kmerLength = 13;
	int kmerSkip = 2;
	int logLevel = 0; //16383, 255
	int minInversionLength = 70; //Advanced
	int minKmerHits = 0;
	int maxCMRs = INT_MAX;
	int maxInitialSegments = 10;
	int maxCLISRuns = 100;
	int readPartCorridor = 40;
	int readPartLength = 256;
	int stdoutMode = 0;
	int subreadaligner = 2;
	int threads = 1;

	float invScoreRatio = 1.0f;
	float scoreMatch = 2.0f;
	float scoreMismatch = -5.0f;
	float scoreGapOpen = -5.0f;
	float scoreGapExtendMax = -5.0f;
	float scoreGapExtendMin = -1.0f;
	float scoreGapDecay = 0.15f;

	bool bam = false;
	bool color = false;
	bool hardClip = false;
	bool log = false;
	bool lowQualitySplit = false;
	bool printAllAlignments = false;
	bool progress = true;
	bool skipSave = false;
	bool smallInversionDetection = true;
	bool updateCheck = false;
	bool verbose = false;
	bool writeUnmapped = true;
	bool nosse = false; //Debug

	char * queryFile = 0;
	char * referenceFile = 0;
	char * outputFile = 0;
	/*
	 * If specified, only reads in the regions specified
	 * by the BED file are read from the input BAM file
	 * (requires BAM input)
	 */
	char * bedFile = 0;
	char * fullCommandLineCall = 0;
	char * vcfFile = 0;


public:

	char const * const getQueryFile() const {
		return queryFile;
	}

	char const * const getReferenceFile() const {
		return referenceFile;
	}

	char const * const getBedFile() const {
		return bedFile;
	}

	char const * const getOutputFile() const {
		return outputFile;
	}

	char const * const getFullCommandLineCall() const {
		return fullCommandLineCall;
	}

	char const * const getVcfFile() const {
		return vcfFile;
	}

	float getMinIdentity() const {
		return minIdentity;
	}

	float getMinResidues() const {
		return minResidues;
	}

	float getSensitivity() const {
		return sensitivity;
	}

	int getBinSize() const {
		return binSize;
	}

	int getCsSearchTableLength() const {
		return csSearchTableLength;
	}

	int getKmerLength() const {
		return kmerLength;
	}

	int getKmerSkip() const {
		return kmerSkip;
	}

	int getLogLevel() const {
		return logLevel;
	}

	int getMinInversionLength() const {
		return minInversionLength;
	}

	int getMinKmerHits() const {
		return minKmerHits;
	}

	int getMaxInitialSegments() const {
		return maxInitialSegments;
	}

	int getMaxCLISRuns() const {
		return maxCLISRuns;
	}

	int getMaxCMRs() const {
		return maxCMRs;
	}

	int getReadPartCorridor() const {
		return readPartCorridor;
	}

	int getReadPartLength() const {
		return readPartLength;
	}

	float getInvScoreRatio() const {
		return invScoreRatio;
	}

	float getScoreMatch() const {
		return scoreMatch;
	}

	float getScoreMismatch() const {
		return scoreMismatch;
	}

	float getScoreGapOpen() const {
		return scoreGapOpen;
	}

	float getScoreExtendMax() const {
		return scoreGapExtendMax;
	}

	float getScoreExtendMin() const {
		return scoreGapExtendMin;
	}

	float getScoreGapDecay() const {
		return scoreGapDecay;
	}

	int getStdoutMode() const {
		return stdoutMode;
	}

	int getSubreadAligner() const {
		return subreadaligner;
	}

	int getThreads() const {
		return threads;
	}

	bool getNoSSE() const {
		return nosse;
	}

	bool getBAM() const {
		return bam;
	}

	bool getColor() const {
		return color;
	}

	bool getHardClip() const {
		return hardClip;
	}

	bool getLog() const {
		return log;
	}

	bool getLowQualitySplit() const {
		return lowQualitySplit;
	}

	bool getSmallInversionDetection() const {
		return smallInversionDetection;
	}

	bool getPrtintAll() const {
		return printAllAlignments;
	}

	bool getProgress() const {
		return progress;
	}

	bool getSkipSave() const {
		return skipSave;
	}

	bool getUpdateCheck() const {
		return updateCheck;
	}

	bool getVerbose() const {
		return verbose;
	}

	bool getWriteUnampped() const {
		return writeUnmapped;
	}

	IConfig() {

	}

	virtual ~IConfig() {
	}

};

#endif //__CONF_H__

typedef void (*pfSetConfig)(IConfig const *);

extern IConfig* _config;
#define Config (*_config)


