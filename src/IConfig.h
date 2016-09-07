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

	int getMinKmerHits() const {
		return minKmerHits;
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

	int getScoreMatch() const {
		return scoreMatch;
	}

	int getScoreMismatch() const {
		return scoreMismatch;
	}

	int getScoreGapOpen() const {
		return scoreGapOpen;
	}

	int getScoreExtend() const {
		return scoreGapExtend;
	}

	int getStdoutMode() const {
		return stdoutMode;
	}

	int getThreads() const {
		return threads;
	}

	bool getFast() const {
		return fast;
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

protected:

	/**
	 * Default values
	 */
	float minIdentity = 0.65f;
	float minResidues = 256.0f;
	float sensitivity = 0.8f;

	int binSize = 4;
	int csSearchTableLength = 0;
	int kmerLength = 13;
	int kmerSkip = 2;
	int logLevel = 0; //16383, 255
	int minKmerHits = 0;
	int maxCMRs = INT_MAX;
	int readPartCorridor = 40;
	int readPartLength = 256;
	int stdoutMode = 0;
	int scoreMatch = 1;
	int scoreMismatch = -4;
	int scoreGapOpen = -2;
	int scoreGapExtend = -2;
	int threads = 1;

	bool bam = false;
	bool color = false;
	bool hardClip = false;
	bool log = false;
	bool lowQualitySplit = true;
	bool printAllAlignments = false;
	bool progress = true;
	bool skipSave = false;
	bool smallInversionDetection = true;
	bool updateCheck = false;
	bool verbose = false;
	bool writeUnmapped = true;
	bool fast = false; //Debug

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

};

#endif //__CONF_H__

typedef void (*pfSetConfig)(IConfig const *);

extern IConfig* _config;
#define Config (*_config)


