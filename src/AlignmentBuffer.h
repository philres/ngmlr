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

#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "GenericReadWriter.h"

#include "ILog.h"
#include "IConfig.h"
#include "intervaltree/IntervalTree.h"
#include "LinearRegression.h"
#include "ConvexAlign.h"
#include "ConvexAlignFast.h"
#include "SAMWriter.h"

#undef module_name
#define module_name "OUTPUT"

//#define TEST_ALIGNER

// Non overlapping part of the reads that are mapped to the reference
// using NextGenMap short-read mapping code (only candidate search and score computation)
struct Anchor {
	// Position of anchor on the read
	int onRead;
	// Position of anchor on the reference
	loc onRef;
	// Alignment score of anchor
	float score;
	// Anchor was mapped to the reverse strand
	bool isReverse;
	// Used for visualization only!
	int type; //0: normal, 1: repetitive, 2: nothing found, 3: socre too low, 4: no coordinates
	// Unique anchors can be used as single intervals, non-unique anchors won't form
	// intervals alone during cLIS
	bool isUnique;
};

static bool sortAnchorOnRead(Anchor a, Anchor b) {
	return a.onRead < b.onRead;
}

static bool sortAnchorOnReadRev(Anchor a, Anchor b) {
	return a.onRead > b.onRead;
}

static bool sortAnchorOnRef(Anchor a, Anchor b) {
	return a.onRef < b.onRef;
}

struct Interval {

public:
	Anchor * anchors;
	int anchorLength;
	int onReadStart;
	int onReadStop;
	loc onRefStart;
	loc onRefStop;
	double m;
	double b;
	double r;
	float score;
	short id;
	bool isReverse;
	bool isProcessed;
	bool isAssigned;

	Interval() {
		anchors = 0;
		anchorLength = 0;
		onReadStart = 0;
		onReadStop = 0;
		onRefStart = 0;
		onRefStop = 0;
		m = 0.0;
		b = 0.0;
		r = 0.0;
		score = 0.0f;
		id = 0;
		isReverse = false;
		isProcessed = false;
		isAssigned = false;

	}

	int lengthOnRead() const {
		return onReadStop - onReadStart;
	}

	loc lengthOnRef() const {
		return llabs(onRefStop - onRefStart);
	}

	virtual ~Interval() {
		if (anchors != 0) {
			delete[] anchors;
			anchors = 0;
			anchorLength = 0;
		}
	}

	void printOneLine() const {
		Log.Message("Interval: %d - %d on read, %lu - %lu on ref, Reverse: %d, Score: %f, Anchors: %d", onReadStart, onReadStop, onRefStart, onRefStop, isReverse, score, anchorLength);

	}

	void print() const {
		Log.Message("Interval with %d anchors on read: %d - %d (diff: %d)", anchorLength, onReadStart, onReadStop, (onReadStop - onReadStart));
		Log.Message("\tAnchor on ref: %lu - %lu (diff: %d)", onRefStart, onRefStop, onRefStop - onRefStart);
		Log.Message("\tReverse: %d, Score: %f", isReverse, score);
	}

private:
	Interval(const Interval & src);

};

class AlignmentBuffer {

private:

	struct Alignment {
		MappedRead * read;
		int scoreId;
	};

	float processTime;
	float alignTime;
//	float overallTime;

	GenericReadWriter* m_Writer;

	static bool first;

	IAlignment * aligner;
#ifdef TEST_ALIGNER
	IAlignment * alignerFast;
#endif

	bool const pacbioDebug;

	bool stdoutPrintDotPlot;
	bool stdoutInversionBed;
	bool stdoutErrorProfile;
	bool printInvCandidateFa;
	bool stdoutPrintMappedSegments;
	bool stdoutPrintAlignCorridor;
	bool stdoutPrintScores;

	IntervalTree::IntervalTree<int> * readCoordsTree;

	int const readPartLength;

	int const maxIntervalNumber;
//	int intervalBufferIndex;
//	Interval ** intervalBuffer;

	void debugAlgnFinished(MappedRead * read);
//	int alignmentCheckForInversion(int const inversionLength,
//			const int refCheckLength, SequenceLocation inversionCheckLocation,
//			uloc inversionMidpointOnRead, const char* const readName,
//			int inversionNumber, char* fullReadSeq);

public:

#define DP_TYPE_UNFILTERED 0
#define DP_TYPE_CLIS 1
#define DP_TYPE_SEQMENTS 200
#define DP_TYPE_SEQMENTS_REG 300
#define DP_TYPE_SEQMENTS_CONS 400
#define DP_TYPE_RESULT 600
#define DP_TYPE_RESULT_CONS 800

#define DP_STATUS_OK 0
#define DP_STATUS_REPETITIVE 1
#define DP_STATUS_NOHIT 2
#define DP_STATUS_LOWSCORE 3
#define DP_STATUS_NOCOORDS 4

#define SV_NONE 0
#define SV_INVERSION 1
#define SV_TRANSLOCATION 2
#define SV_UNKNOWN 3

	// A list of intervals that are "compatible" meaning they are located in
	// a "corridor" that is small enough to be handled by the alignment
	// algorithm
	// TODO: remove fixed length!
	struct MappedSegment {

		static int const maxLength = 1000;

		Interval * list[maxLength];
		size_t length;

		MappedSegment() {
			length = 0;
		}
	};

	// Signed position on the reference
	// For cLIS reverse positions are represented as negative
	// numbers (TODO: check whether this in neccessary!)
	typedef long long loc;

	int * cLIS(Anchor * anchors, int const anchorsLenght, int & lisLength);

	void addAnchorAsInterval(Anchor const & anchor, MappedSegment & segment);
	Interval * toInterval(Anchor const & anchor);

	bool isCompatible(Interval const * a, Interval const * b);
	bool isCompatible(Anchor const & anchor, MappedSegment const & segment);
	bool isCompatible(Anchor const & anchor, Interval const * interval);
	bool isContained(Interval const * a, Interval const * b);
	bool isSameDirection(Interval const * a, Interval const * b);
	bool isDuplication(Interval const *, int, Interval const *);

	/**
	 * Returns the number of overlapping read bps between the two intervals
	 */
	int getOverlapOnRead(Interval const * a, Interval const * b);

	/**
	 * Prints message if --verbose is specified.
	 * Adds <tabs> tabstops before printing
	 */
	void verbose(int const tabs, bool const newLine, char const * const s, ...);

	/**
	 * Print interval if --verbose is specified
	 */
	void verbose(int const tabs, char const * const s, Interval * interval);

	/**
	 * Removes <readBp> bp from the start of an interval.
	 * RefBp are estimated from the ratio of total ref and read bp in
	 * the interval.
	 */
	bool shortenIntervalStart(Interval * interval, int const readBp, bool readOnly = false);

	/**
	 * Removes <readBp> bp from the end of an interval.
	 * RefBp are estimated from the ratio of total ref and read bp in
	 * the interval.
	 */
	bool shortenIntervalEnd(Interval * interval, int const readBp, bool readOnly = false);


	Interval * mergeIntervals(Interval * a, Interval * b);
	Interval * * infereCMRsfromAnchors(int & intervalsIndex,
			Anchor * allFwdAnchors, int allFwdAnchorsLength,
			Anchor * allRevAnchors, int allRevAnchorsLength, MappedRead * read);
	Interval * * getIntervalsFromAnchors(int & intervalsIndex,
			Anchor * allFwdAnchors, int allFwdAnchorsLength, Anchor * allRevAnchors,
			int allRevAnchorsLength, MappedRead * read);


	Align computeAlignment(MappedRead* read, int const scoreId,
			int const corridor);

	Align * computeAlignment(Interval const * interval, int const corridor,
			char * const readSeq, size_t const readLength, int const QStart,
			int const QEnd, int fullReadLength, MappedRead const * const read,
			bool const realign, bool const fullAlignment, bool const shortRead);

	int estimateCorridor(Interval const * interval);
	char * const extractReadSeq(int const readSeqLen,
			Interval const * interval, MappedRead* read, bool const revComp);

	Align * alignInterval(MappedRead const * const read,
			Interval const * interval, char * const readSeq,
			size_t const readSeqLen, bool const realign, bool const fullAlignment);
	void alignSingleOrMultipleIntervals(MappedRead * read,
			Interval const * const interval, LocationScore * tmp,
			Align * tmpAling, int & alignIndex);

	int realign(int const svType, Interval const * interval,
			Interval * leftOfInv, Interval * rightOfInv, MappedRead * read,
			Align * tmpAling, int & alignIndex, LocationScore * tmp, int mq);

	bool constructMappedSegements(MappedSegment * segments, Interval * interval,
			size_t & segmentsIndex);

	bool reconcileRead(ReadGroup * group);

	Interval * * consolidateSegments(MappedSegment * segments,
			size_t segmentsIndex, int & intervalsIndex);
	void consolidateSegment(Interval * interval, int & intervalsIndex,
			MappedSegment segment);

	int detectMisalignment(Align const * const align, Interval const * interval,
			char * readPartSeq, Interval * leftOfInv, Interval * rightOfInv,
			MappedRead * read);

	int checkForSV(Align const * const align, Interval const * interval, char * fullReadSeq, uloc inversionMidpointOnRef, uloc inversionMidpointOnRead, int inversionLength, MappedRead * read);

	void resolveReadOverlap(Interval * first, Interval * second);

	void processLongRead(ReadGroup * group);
	void processLongReadLIS(ReadGroup * group);
	void processShortRead(MappedRead * read);

	int computeMappingQuality(Align const & alignment, int readLength);

	void printDotPlotLine(int const id, char const * const name,
			int const onReadStart, int const onReadStop, loc const onRefStart,
			loc const onRefStop, float const score, bool const isReverse,
			int const type, int const status);

	void printDotPlotLine(int const id, char const * const name,
	REAL const m, REAL const b, REAL const r, float const score,
			bool const isReverse, int const type, int const status);

	AlignmentBuffer(const char* const filename) :
			pacbioDebug(Config.getVerbose()), readCoordsTree(0), readPartLength(Config.getReadPartLength()), maxIntervalNumber(Config.getMaxInitialSegments()) {

		m_Writer = 0;

		m_Writer = (GenericReadWriter*) new SAMWriter((FileWriter*) NGM.getWriter());

		if (first) {
			m_Writer->WriteProlog();
			first = false;
		}

		processTime = 0.0f;
//		overallTime = 0.0f;
		alignTime = 0.0f;

		if (first) {
			m_Writer->WriteProlog();
			first = false;
		}

		stdoutPrintDotPlot = Config.getStdoutMode() == 1;
		stdoutInversionBed = Config.getStdoutMode() == 2;
		stdoutErrorProfile = Config.getStdoutMode() == 3;
		printInvCandidateFa = Config.getStdoutMode() == 4;
		stdoutPrintMappedSegments = Config.getStdoutMode() == 5;
		stdoutPrintAlignCorridor = Config.getStdoutMode() == 6;
		stdoutPrintScores = Config.getStdoutMode() == 7;

		Log.Verbose("Alignment batchsize = %i", batchSize);

		//	IAlignment * aligner = new StrippedSW();
		if (Config.getNoSSE()) {
			aligner = new Convex::ConvexAlign(
								Config.getStdoutMode(),
								Config.getScoreMatch(),
								Config.getScoreMismatch(),
								Config.getScoreGapOpen(),
								Config.getScoreExtendMax(),
								Config.getScoreExtendMin(),
								Config.getScoreGapDecay());
		} else {
			aligner = new Convex::ConvexAlignFast(
								Config.getStdoutMode(),
								Config.getScoreMatch(),
								Config.getScoreMismatch(),
								Config.getScoreGapOpen(),
								Config.getScoreExtendMax(),
								Config.getScoreExtendMin(),
								Config.getScoreGapDecay());
		}
#ifdef TEST_ALIGNER
		alignerFast = new Convex::ConvexAlignFast(0);
#endif
	}

	virtual ~AlignmentBuffer() {
		delete m_Writer;
		m_Writer = 0;

//		delete[] intervalBuffer;
//		intervalBuffer = 0;

		delete aligner;
		aligner = 0;

#ifdef TEST_ALIGNER
		delete alignerFast;
		alignerFast = 0;
#endif
	}

	void DoRun();

	int GetStage() const {
		return 4;
	}

	inline const char* GetName() const {
		return "Output";
	}

	void addRead(MappedRead * read, int scoreID);
	void flush();

	float getTime() {
//		float tmp = overallTime;
//		overallTime = 0;
		return 0;
	}

	float getProcessTime() {
		float tmp = processTime;
//		processTime = 0;
		return tmp;
	}

	float getAlignTime() {
		float tmp = alignTime;
//		alignTime = 0;
		return tmp;
	}

	void SaveRead(MappedRead* read, bool mapped = true);
	void WriteRead(MappedRead* read, bool mapped);
}
;

#endif
