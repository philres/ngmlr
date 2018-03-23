/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __ICONF_H__
#define __ICONF_H__

#include "Types.h"

#include <climits>
#include <cstring>

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

	int maxSegmentNumberPerKb = 1;

	int maxCLISRuns = 100;

	int readPartCorridor = 40;
	int readPartLength = 256;
	int stdoutMode = 0;
	int subreadaligner = 2;
	int threads = 1;
	int maxReadNameLength = 500;

	ulong maxMatrixSizeMB = 10000;

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
	bool bamCigarFix = false;
	bool skipAlign = false;

	char * queryFile = 0;
	char * referenceFile = 0;
	char * outputFile = 0;
	char * rgId = 0;
	char * rgSm = 0;
	char * rgLb = 0;
	char * rgPl = 0;
	char * rgDs = 0;
	char * rgDt = 0;
	char * rgPu = 0;
	char * rgPi = 0;
	char * rgPg = 0;
	char * rgCn = 0;
	char * rgFo = 0;
	char * rgKs = 0;

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

	char const * const getRgId() const {
		return rgId;
	}

	char const * const getRgSm() const {
		return rgSm;
	}

	char const * const getRgLb() const {
		return rgLb;
	}

	char const * const getRgPl() const {
		return rgPl;
	}

	char const * const getRgDs() const {
		return rgDs;
	}

	char const * const getRgDt() const {
		return rgDt;
	}

	char const * const getRgPu() const {
		return rgPu;
	}

	char const * const getRgPi() const {
		return rgPi;
	}

	char const * const getRgPg() const {
		return rgPg;
	}

	char const * const getRgCn() const {
		return rgCn;
	}

	char const * const getRgFo() const {
		return rgFo;
	}

	char const * const getRgKs() const {
		return rgKs;
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

	int getMaxSegmentNumberPerKb(int const readLength) const {
		//int maxSplits = int((read->length / 1000.0) * Config.getMaxSegmentNumberPerKb() + 0.5);
		int maxSegments = (int)((readLength / 1000.0) * maxSegmentNumberPerKb + 0.5);
		return maxSegments < 1 ? 1 : maxSegments;
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

	int getMaxReadNameLength() const {
		return maxReadNameLength;
	}

	ulong getMaxMatrixSizeMB() const {
		return maxMatrixSizeMB;
	}

	bool getNoSSE() const {
		return nosse;
	}

	bool getBamCigarFix() const {
		return bamCigarFix;
	}

	bool getSkipalign() const {
		return skipAlign;
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

	bool isStdIn() const {
		return strcmp(getQueryFile(), "/dev/stdin") == 0;
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


