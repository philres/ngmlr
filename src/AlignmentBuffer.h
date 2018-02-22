/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "GenericReadWriter.h"

#include <memory>

#include "ILog.h"
#include "IConfig.h"
#include "intervaltree/IntervalTree.h"
#include "LinearRegression.h"
#include "ConvexAlign.h"
#include "ConvexAlignFast.h"
#include "StrippedSW.h"
#include "SAMWriter.h"

#undef module_name
#define module_name "OUTPUT"


using std::unique_ptr;
//#define TEST_ALIGNER

static bool sortAnchorOnRead(Anchor a, Anchor b) {
	return a.onRead < b.onRead;
}

static bool sortAnchorOnReadRev(Anchor a, Anchor b) {
	return a.onRead > b.onRead;
}

static bool sortAnchorOnRef(Anchor a, Anchor b) {
	return a.onRef < b.onRef;
}


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
	IAlignment * overlapCheckAligner;

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

	int maxSegmentCount;
//	float const maxIntervalNumberPerKb;
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

	int * cLIS(Anchor * anchors, int const anchorsLenght, int & lisLength);

	/**
	 * Checks if shorter interval is long enough to span deletion/insertion
	 * (if present). Only if intervals are close on read or on reference.
	 * Returns true if there is a large gap between intervals
	 */
	bool canSpanDeletionInsertion(Interval const * a, Interval const * b, REAL corridorSize);

	/**
	 * Checks whether the gap between two intervals spans a chromosome border
	 */
	bool spansChromosomeBorder(Interval const * a, Interval const * b);

	void addAnchorAsInterval(Anchor const & anchor, MappedSegment & segment);
	Interval * toInterval(Anchor const & anchor);

	/**
	 * Test if interval is contained in corridor
	 */
	bool isIntervalInCorridor(REAL k, REAL d, REAL corridor,
			Interval const * testee, bool const switched);

	bool isCompatible(Interval const * a, Interval const * b, REAL corridorSize = 8192.0);
//	bool isCompatible(Anchor const & anchor, MappedSegment const & segment);
//	bool isCompatible(Anchor const & anchor, Interval const * interval);
	bool isContained(Interval const * a, Interval const * b);
	bool isSameDirection(Interval const * a, Interval const * b);
	bool isDuplication(Interval const *, Interval const *, loc & dupLength);

	/**
	 * Distance between two intervals on read
	 */
	int getDistanceOnRead(Interval const * a, Interval const * b);
	/**
	 * Distance between two intervals on ref
	 */
	loc getDistanceOnRef(Interval const * a, Interval const * b);
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
	void verbose(int const tabs, char const * const s, Interval const * const interval);

	/**
	 * Removes <readBp> bp from the start of an interval.
	 * RefBp are estimated from the ratio of total ref and read bp in
	 * the interval.
	 */
	bool shortenIntervalStart(Interval * interval, int const readBp);

	/**
	 * Removes <readBp> bp from the end of an interval.
	 * RefBp are estimated from the ratio of total ref and read bp in
	 * the interval.
	 */
	bool shortenIntervalEnd(Interval * interval, int const readBp);

	bool extendIntervalStart(Interval * interval, int const readBp);

	bool extendIntervalStop(Interval * interval, int const readBp, int const readLength);

	/**
	 * Check if gap between first and second interval overlaps with another interval
	 */
	bool gapOverlapsWithInterval(Interval * first, Interval * second, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read);

	/**
	 * Checks if interval overlaps with any other interval in interval tree
	 */
	bool gapOverlapsWithInterval(Interval * gap, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read);

	/**
	 * Checks if gap between second and read end overlaps with any other interval
	 */
	bool gapToEndOverlapsWithInterval(Interval * second, int const readLength, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read);

	/**
	 * Checks if gap between read start and second overlaps with any other interval
	 */
	bool gapFromStartOverlapsWithInterval(Interval * second, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read);

	/**
	 * Extends both intervals to close the gap on the read
	 */
	void closeGapOnRead(Interval * first, Interval * second, int const readLength);

	/**
	 * Extends interval to read start, if not overlapping with other interval
	 */
	void extendToReadStart(Interval * interval, int const readLength, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read);

	/**
	 * Extends interval to read stop, if not overlapping with other interval
	 */
	void extendToReadStop(Interval * interval, int const readLength, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read);


	/**
	 * Compute alignment band with fixed length from anchors
	 */
	CorridorLine * getCorridorEndpointsWithAnchors(Interval const * interval, int const corridorMultiplier, char const * refSeq, char const * readSeq, int & corridorHeight, int const externalQStart, int const readPartLength,
			int const fullReadLength, bool const realign);

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
			char const * const readSeq, size_t const readLength, int const QStart,
			int const QEnd, int fullReadLength, MappedRead const * const read,
			bool const realign, bool const fullAlignment, bool const shortRead);

	int estimateCorridor(Interval const * interval);

	unique_ptr<char []> extractReadSeq(int const readSeqLen,
			Interval const * interval, MappedRead* read, bool const revComp = false);

	unique_ptr<char []> extractReadSeq(int const readSeqLen,
			int const onReadStart, bool const isReverse, MappedRead* read,
			bool const revComp);

	Align * alignInterval(MappedRead const * const read,
			Interval const * interval, char const * const readSeq,
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

	Interval * getIntervalFromAlign(Align const * const align, LocationScore const * const score, int const i, int const readLength);

	Interval * * consolidateSegments(MappedSegment * segments,
			size_t segmentsIndex, int & intervalsIndex);
	void consolidateSegment(Interval * interval, int & intervalsIndex,
			MappedSegment segment);

	int detectMisalignment(Align const * const align, Interval const * interval,
			char const * const readPartSeq, Interval * leftOfInv, Interval * rightOfInv,
			MappedRead * read);

	int checkForSV(Align const * const align, Interval const * interval, char const * const fullReadSeq, uloc inversionMidpointOnRef, uloc inversionMidpointOnRead, int inversionLength, MappedRead * read);

	/**
	 * Align interval using StrippedSW and return score
	 */
	float scoreInterval(Interval * interval, MappedRead * read);

	/**
	 * Extracts sequence from reference genome
	 */
	char const * const extractReferenceSequenceForAlignment(Interval const*& interval, int & refSeqLength);

	/**
	 * Extracts sequence from reference genome
	 */
	char const * const extractReferenceSequenceForAlignment(loc const onRefStart, loc const onRefStop, int & refSeqLength);

	void processLongReadLIS(ReadGroup * group);
	void processShortRead(MappedRead * read);

	int computeMappingQuality(Align const & alignment, int readLength);

	void printDotPlotLine(int const id, char const * const name, int const onReadStart, int const onReadStop, loc const onRefStart, loc const onRefStop, float const score, bool const isReverse, int const type, int const status);

	void printDotPlotLine(int const id, char const * const name,
	REAL const m, REAL const b, REAL const r, float const score, bool const isReverse, int const type, int const status);

	AlignmentBuffer(const char* const filename) :
			pacbioDebug(Config.getVerbose()), readCoordsTree(0), readPartLength(Config.getReadPartLength()), maxSegmentCount(0)/*, maxIntervalNumberPerKb(Config.getMaxSegmentNumberPerKb())*/ {

		m_Writer = (GenericReadWriter*) new SAMWriter((FileWriter*) NGM.getWriter());

		if (first) {
			m_Writer->WriteProlog();
			first = false;
		}

		processTime = 0.0f;
//		overallTime = 0.0f;
		alignTime = 0.0f;

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

		overlapCheckAligner = new StrippedSW();
	}

	virtual ~AlignmentBuffer() {
		delete m_Writer;
		m_Writer = nullptr;

		delete aligner;
		aligner = nullptr;

#ifdef TEST_ALIGNER
		delete alignerFast;
		alignerFast = nullptr;
#endif
		if(overlapCheckAligner != nullptr) {
			delete overlapCheckAligner;
			overlapCheckAligner = nullptr;
		}
	}

	int GetStage() const {
		return 4;
	}

	inline const char* GetName() const {
		return "Output";
	}

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
