/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __IALIGNMENT_H__
#define __IALIGNMENT_H__

#include <stdlib.h>

// Signed position on the reference
// For cLIS reverse positions are represented as negative
// numbers (TODO: check whether this in neccessary!)
typedef long long loc;


struct PositionNM {

	PositionNM() :
			refPosition(0), readPosition(0), nm(0) {

	}

	int refPosition;
	int readPosition;
	int nm;

};

struct CorridorLine {
	int offset;
	int length;
	unsigned long offsetInMatrix;
};

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

private:
	Interval(const Interval & src);

};



struct Align {

public:

//	static volatile int sInstanceCount;

	Align() :
			pBuffer1(0), pBuffer2(0), mappedInterval(0), nmPerPosition(0), nmPerPostionLength(0), alignmentLength(
					0), PositionOffset(0), QStart(0), QEnd(0), Score(0.0f), Identity(
					0.0f), NM(0), MQ(0), cigarOpCount(0), maxBufferLength(20000), maxMdBufferLength(20000), skip(false), primary(false), svType(0) {

	}

	virtual ~Align() {

	}

	char * pBuffer1; // = pCigar = pRef
	char * pBuffer2; // = pMD = pQry
	PositionNM * nmPerPosition;
//	void * ExtendedData;
	Interval * mappedInterval;
	PositionNM firstPosition;
	PositionNM lastPosition;
	int nmPerPostionLength;
	int alignmentLength;
	int PositionOffset; // Position in Ref, an der das Alignment beginnt
	int QStart; // Anzahl Basen, die beim Qry am Anfang abgeschnitten wurden
	int QEnd; // Anzahl Basen, die beim Qry am Ende abgeschnitten wurden
	float Score;
	float Identity;
	int NM;
	int MQ;
	int cigarOpCount;
	int maxBufferLength;
	int maxMdBufferLength;
	bool skip;
	bool primary;
	/**
	 * Bitflag:
	 * 		0x1 - Start clipping because of N
	 * 		0x2 - End clipping because of N
	 */
	int svType;

	void setBitFlag(int const i) {
		svType = svType | i;
	}

	int getAlignedReadBp(int const readLength) {
		return readLength - QStart - QEnd;
	}

	void clearBuffer() {
		if (pBuffer1 != 0) {
			delete[] pBuffer1;
			pBuffer1 = 0;
		}
		if (pBuffer2 != 0) {
			delete[] pBuffer2;
			pBuffer2 = 0;
		}
		if (mappedInterval != 0) {
			delete[] mappedInterval;
			mappedInterval = 0;
		}
	}

	void clearNmPerPosition() {
		if (nmPerPosition != 0) {
			delete[] nmPerPosition;
			nmPerPosition = 0;
			nmPerPostionLength = 0;
		}
	}

//private:
//	Align(const Align & src);

};

static int const cCookie = 0x10201130;

/*
 Anmerkung zum Parameter mode:

 int AlignmentType = mode & 0xFF;	// 0..Smith-Waterman, 1..Needleman-Wunsch
 int ReportType = (mode >> 8) & 0xFF;	// 0..Plain alignment (Ref+Qry), 1..SAM (Cigar+MD)
 bool BSMappingActive = mode & 0x10000;

 Anmerkung BS-Mapping:

 extData zeigt bei BSMappingActive == true auf ein Flag-Array (char*) der L�nge batchSize,
 wobei bei 0 die TC-Match-Funktion, bei 1 die AG-Match-Funktion verwendet werden soll:

 if (extData[i] == 0) -> TC-Matching f�r ref/qry-Paar i
 if (extData[i] == 1) -> AG-Matching     - "" -

 */
class IAlignment {
public:
	virtual int GetScoreBatchSize() const = 0;
	virtual int GetAlignBatchSize() const = 0;

	virtual int BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, float * const results,
			void * extData) = 0;

	virtual int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData) {
		return 0;
	}

	virtual int SingleAlign(int const mode, CorridorLine * corridor,
			int const corridorHeight, char const * const refSeq,
			char const * const qrySeq, Align & result, int const externalQStart,
			int const externalQEnd, void * extData) {
		return 0;
	}

	virtual int SingleScore(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			float & result, void * extData) {
		return 0;
	}

	virtual int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results,
			void * extData) = 0;

	virtual ~IAlignment() {
	}
};

typedef IAlignment * (*pfCreateAlignment)(int const gpu_id);
typedef void (*pfDeleteAlignment)(IAlignment*);

#endif
