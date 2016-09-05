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
		/**
		 * Default values
		 */
		minIdentity = 0.65f;
		minResidues = 50.0f;
		sensitivity = 0.8f;


		binSize = 4;
		csSearchTableLength = 0;
		kmerLength = 13;
		kmerSkip = 2;
		logLevel = 0; //16383, 255
		minKmerHits = 0;
		maxCMRs = INT_MAX;
		readPartCorridor = 40;
		readPartLength = 256;
		stdoutMode = 0;
		scoreMatch = 1;
		scoreMismatch = -4;
		scoreGapOpen = -2;
		scoreGapExtend = -2;
		threads = 1;


		bam = false;
		color = false;
		hardClip = false;
		log = false;
		lowQualitySplit = true;
		printAllAlignments = false;
		progress = true;
		skipSave = false;
		smallInversionDetection = true;
		updateCheck = false;
		verbose = false;
		writeUnmapped = true;
		fast = false; //Debug


		queryFile = 0;
		referenceFile = 0;
		outputFile = 0;
		bedFile = 0;
		fullCommandLineCall = 0;
		vcfFile = 0;
	}

	virtual ~IConfig() {
	}

protected:

	/**
	 * Float parameters
	 */

	float minIdentity;
	float minResidues;
	float sensitivity;

	/**
	 * Int parameters
	 */

	int binSize;

	int csSearchTableLength;

	int kmerLength;

	int kmerSkip;

	int minKmerHits;

	int logLevel;

	int maxCMRs;

	int readPartCorridor;

	int readPartLength;

	int scoreMatch;

	int scoreMismatch;

	int scoreGapOpen;

	int scoreGapExtend;

	int stdoutMode;

	int threads;

	/**
	 * Bool parameters
	 */
	bool bam;

	bool color;

	bool fast; //Debug

	bool hardClip;

	bool log;

	bool lowQualitySplit;

	bool smallInversionDetection;

	bool printAllAlignments; //Debug

	bool progress;

	bool skipSave;

	bool updateCheck;

	bool verbose;

	bool writeUnmapped;

	/**
	 * String parameters
	 */

	/*
	 * If specified, only reads in the regions specified
	 * by the BED file are read from the input BAM file
	 * (requires BAM input)
	 */
	char * bedFile;

	char * fullCommandLineCall;

	char * outputFile;

	char * queryFile;

	char * referenceFile;

	char * vcfFile;



};

#endif //__CONF_H__

typedef void (*pfSetConfig)(IConfig const *);

extern IConfig* _config;
#define Config (*_config)


