#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "GenericReadWriter.h"

#include "SAMWriter.h"
#include "BAMWriter.h"
#include "ScoreWriter.h"
#include "StrippedSW.h"

#include "intervaltree/IntervalTree.h"

#undef module_name
#define module_name "OUTPUT"

class AlignmentBuffer {

private:

	struct Alignment {
		MappedRead * read;
		int scoreId;
	};

	int const outputformat;
	int const alignmode;
	bool m_EnableBS;
	int const batchSize;
	int const corridor;
	uloc const refMaxLen;
	int const min_mq;

	Alignment * reads;
	int nReads;
	char const * * qryBuffer;
	char const * * refBuffer;
	char * m_DirBuffer;
	int dbLen;
	Align * alignBuffer;
	char * dBuffer;
	char * dummy;

	long pairInsertCount;
	long pairInsertSum;
	long brokenPairs;

	float processTime;
	float alignTime;
	float overallTime;

	GenericReadWriter* m_Writer;

	static bool first;

	IAlignment * aligner;

	bool const argos;

	bool const pacbioDebug;

	IntervalTree::IntervalTree<int> * readCoordsTree;

	void debugAlgnFinished(MappedRead * read);

public:

#define TYPE_UNFILTERED 0
#define TYPE_LIS_SEQAN 1
#define TYPE_CLIS 200
#define TYPE_RESULT 1000

#define STATUS_OK 0
#define STATUS_REPETITIVE 1
#define STATUS_NOHIT 2
#define STATUS_LOWSCORE 3

	struct Interval {
		int onReadStart;
		int onReadStop;
		loc onRefStart;
		loc onRefStop;
		float score;
		bool isReverse;

		void printOneLine() {
			Log.Message("Interval: %d - %d on read, %lu - %lu on ref, Reverse: %d, Score: %f", onReadStart, onReadStop, onRefStart, onRefStop, isReverse, score);
//			Log.Message("\tAnchor on ref: %lu - %lu (diff: %d)", onRefStart, onRefStop, onRefStop - onRefStart);
//			Log.Message("\tReverse: %d, Score: %f", isReverse, score);
		}

		void print() {
			Log.Message("Interval on read: %d - %d (diff: %d)", onReadStart, onReadStop, (onReadStop - onReadStart));
			Log.Message("\tAnchor on ref: %lu - %lu (diff: %d)", onRefStart, onRefStop, onRefStop - onRefStart);
			Log.Message("\tReverse: %d, Score: %f", isReverse, score);
		}

		void print(int id, char * name, int intNumber) {
//			printf("%d\t%s\t%d\t%d\t%llu\t%llu\t%f\t%d\t%d\t%d\n", id, name, onReadStart, onReadStop, onRefStart, onRefStop,
//					score, isReverse, TYPE_RESULT + intNumber, STATUS_OK);
		}

	};

	struct MappedSegment {
		Interval list[100];
		size_t length;
	};

	typedef long long loc;

	struct Anchor {
		int onRead;
		loc onRef;
		float score;
		bool isReverse;
		int type; //0: normal, 1: repetitive, 2: nothing found, 3: socre too low
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

	void getSeqAnLIS(int allFwdAnchorsLength, int id,
			Anchor* allFwdAnchors, char* name);

	//Type: 0 unfiltered, 1 filtered
	void printDotPlot(char * name, int id, Anchor * anchors, int anchorsLength,
			int filtered) {
		for (int i = 0; i < anchorsLength; ++i) {
//			printf("%d\t%s\t%d\t%d\t%llu\t%llu\t%f\t%d\t%d\t%d\n", id, name,
//					anchors[i].onRead * 512, anchors[i].onRead * 512 + 512,
//					anchors[i].onRef, anchors[i].onRef + 512, anchors[i].score,
//					anchors[i].isReverse, filtered, anchors[i].type);
		}
	}

	static ulong alignmentCount;

	int * cLIS(Anchor * anchors, int const anchorsLenght, int & lisLength);

	bool isCompatible(Interval a, Interval b);
	bool isContained(Interval a, Interval b);
	bool isContainedOnRead(Interval a, Interval b);
	bool isSameDirection(Interval a, Interval b);
	Interval splitInterval(Interval a,
			Interval b);
	Interval mergeIntervals(Interval a, Interval b);
	Interval * findSubsequences(char * name, int id, Anchor * allFwdAnchors,
			int allFwdAnchorsLength, Anchor * allRevAnchors,
			int allRevAnchorsLength, int readParts, int readLenth,
			int & intervalsIndex);

	Align computeAlignment(MappedRead* read, int const scoreId,
			int const corridor);

	Align computeAlignment(uloc const position, int const corridor,
			char * const readSeq, size_t const readLength, int const QStart, int const QEnd, int fullReadLength);

	int estimateCorridor(const Interval & interval);
	char * const extractReadSeq(const size_t& readSeqLen, Interval & interval,
			MappedRead* read);

	Align alignInterval(MappedRead const * const read, Interval const interval,
			char * const readSeq, size_t const readSeqLen);
	void alignSingleOrMultipleIntervals(MappedRead * read, Interval interval, LocationScore * tmp, Align * tmpAling, int & alignIndex);

	bool alignInversion(Interval interval, Interval leftOfInv, Interval inv,
			Interval rightOfInv, MappedRead * read, Align * tmpAling,
			int & alignIndex, LocationScore * tmp, int mq);

//	bool constructMappedSegements(Interval * intervals,
//			Interval interval, int & intervalsIndex);

	bool constructMappedSegements(MappedSegment * segments,
			Interval interval, size_t & segmentsIndex);

//	bool sortIntervalsInSegment(Interval a, Interval b);

	void reconcileRead(ReadGroup * group);

	Interval * consolidateSegments(MappedSegment * segments, size_t segmentsIndex, int & intervalsIndex);
	void consolidateSegment(Interval * interval, int & intervalsIndex, MappedSegment segment);

	bool inversionDetection(Align const align, Interval const interval, int const length,
			char * fullReadSeq, Interval & leftOfInv, Interval & rightOfInv, Interval & inv, char const * const readName);

	void processLongRead(ReadGroup * group);
	void processLongReadLIS(ReadGroup * group);

	int computeMappingQuality(Align const & alignment, int readLength);

	AlignmentBuffer(const char* const filename, IAlignment * mAligner) :
	batchSize(mAligner->GetAlignBatchSize() / 2), outputformat(
			NGM.GetOutputFormat()),
	alignmode(Config.GetInt(MODE, 0, 1)),
	corridor(Config.GetInt("corridor")),
	refMaxLen((Config.GetInt("qry_max_len") + corridor) | 1 + 1),
	min_mq(Config.GetInt(MIN_MQ)),
	aligner(mAligner), argos(Config.Exists(ARGOS)), pacbioDebug(Config.GetInt(PACBIOLOG) == 1), readCoordsTree(0) {
		pairInsertSum = 0;
		pairInsertCount = 0;
		brokenPairs = 0;
		m_Writer = 0;
		nReads = 0;

		m_EnableBS = false;
		m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

		int const outputformat = NGM.GetOutputFormat();

		if(Config.Exists(ARGOS)) {
			m_Writer = (GenericReadWriter*) new ScoreWriter((FileWriter*)NGM.getWriter());
		} else {
			switch (outputformat) {
				case 0:
				Log.Error("This output format is not supported any more.");
				Fatal();
				break;
				case 1:
				m_Writer = (GenericReadWriter*) new SAMWriter((FileWriter*)NGM.getWriter());
				break;
				case 2:
				m_Writer = (GenericReadWriter*) new BAMWriter((FileWriterBam*)NGM.getWriter(), filename);
				break;
				default:
				break;
			}
		}

		if(first) {
			m_Writer->WriteProlog();
			first = false;
		}

		reads = new Alignment[batchSize];

		qryBuffer = new char const *[batchSize];
		refBuffer = new char const *[batchSize];

		for (int i = 0; i < batchSize; ++i) {
			refBuffer[i] = new char[refMaxLen];
		}

		m_DirBuffer = new char[batchSize];

		alignBuffer = new Align[batchSize];
		dbLen = std::max(1, Config.GetInt("qry_max_len")) * 8;
		dBuffer = new char[dbLen];

		dummy = new char[refMaxLen];
		memset(dummy, '\0', refMaxLen);
		//dummy[Config.GetInt("qry_max_len") - 1] = '\0';

		processTime = 0.0f;
		overallTime = 0.0f;
		alignTime = 0.0f;
		//}
		//}

		if(first) {
			m_Writer->WriteProlog();
			first = false;
		}

		Log.Verbose("Alignment batchsize = %i", batchSize);

	}

	virtual ~AlignmentBuffer() {
		delete m_Writer;
		delete[] m_DirBuffer;
		m_DirBuffer = 0;

		delete[] dummy;
		dummy = 0;

		for (int i = 0; i < batchSize; ++i) {
			delete[] refBuffer[i];
			refBuffer[i] = 0;
		}
		delete[] qryBuffer;
		delete[] refBuffer;
		delete[] alignBuffer;

		delete[] reads;
		delete[] dBuffer;

		//m_Writer->WriteEpilog();

		//delete m_Writer;
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
		float tmp = overallTime;
		overallTime = 0;
		return tmp;
	}

	float getProcessTime() {
		float tmp = processTime;
		processTime = 0;
		return tmp;
	}

	float getAlignTime() {
		float tmp = alignTime;
		alignTime = 0;
		return tmp;
	}

	void SaveRead(MappedRead* read, bool mapped = true);
	void WriteRead(MappedRead* read, bool mapped);
};

#endif
