/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
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

	float scoreMatch = 2.0f;
	float scoreMismatch = -5.0f;
	float scoreGapOpen = -5.0f;
	float scoreGapExtendMax = -5.0f;
	float scoreGapExtendMin = -1.0f;
	float scoreGapDecay = 0.15f;
	float threads = 1;

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

	auto getMinIdentity() const {
		return minIdentity;
	}

	auto getMinResidues() const {
		return minResidues;
	}

	auto getSensitivity() const {
		return sensitivity;
	}

	auto getBinSize() const {
		return binSize;
	}

	auto getCsSearchTableLength() const {
		return csSearchTableLength;
	}

	auto getKmerLength() const {
		return kmerLength;
	}

	auto getKmerSkip() const {
		return kmerSkip;
	}

	auto getLogLevel() const {
		return logLevel;
	}

	auto getMinInversionLength() const {
		return minInversionLength;
	}

	auto getMinKmerHits() const {
		return minKmerHits;
	}

	auto getMaxInitialSegments() const {
		return maxInitialSegments;
	}

	auto getMaxCLISRuns() const {
		return maxCLISRuns;
	}

	auto getMaxCMRs() const {
		return maxCMRs;
	}

	auto getReadPartCorridor() const {
		return readPartCorridor;
	}

	auto getReadPartLength() const {
		return readPartLength;
	}

	auto getScoreMatch() const {
		return scoreMatch;
	}

	auto getScoreMismatch() const {
		return scoreMismatch;
	}

	auto getScoreGapOpen() const {
		return scoreGapOpen;
	}

	auto getScoreExtendMax() const {
		return scoreGapExtendMax;
	}

	auto getScoreExtendMin() const {
		return scoreGapExtendMin;
	}

	auto getScoreGapDecay() const {
		return scoreGapDecay;
	}

	auto getStdoutMode() const {
		return stdoutMode;
	}

	auto getSubreadAligner() const {
		return subreadaligner;
	}

	auto getThreads() const {
		return threads;
	}

	auto getNoSSE() const {
		return nosse;
	}

	auto getBAM() const {
		return bam;
	}

	auto getColor() const {
		return color;
	}

	auto getHardClip() const {
		return hardClip;
	}

	auto getLog() const {
		return log;
	}

	auto getLowQualitySplit() const {
		return lowQualitySplit;
	}

	auto getSmallInversionDetection() const {
		return smallInversionDetection;
	}

	auto getPrtintAll() const {
		return printAllAlignments;
	}

	auto getProgress() const {
		return progress;
	}

	auto getSkipSave() const {
		return skipSave;
	}

	auto getUpdateCheck() const {
		return updateCheck;
	}

	auto getVerbose() const {
		return verbose;
	}

	auto getWriteUnampped() const {
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


