#ifndef __MAPPEDREAD_H__
#define __MAPPEDREAD_H__

#include "LocationScore.h"

#include <memory.h>
#include <vector>
#include <algorithm>
#include "Types.h"
#include "SequenceProvider.h"
#include "ILog.h"
#include "ReadStatus.h"

using NGMNames::ReadStatus;
// A read with a list of scores fitting the initial predicate
// (i.e., each score is greater than cs_threshhold)
struct MappedRead {
private:

public:
	int const ReadId;

	//Calculated is initialized with -1. This shows that the read has not passed the candidate search
	//When the read is submitted to the score computation calculated is set to 0. Each time a score
	//is computed for this read, Calculated is increased by one
	volatile int Calculated; // Amount of Scores updated by SW

	int const qryMaxLen;

	int TopScore; // Position of TopScore
	volatile int Lock; // sync for paired read selection
	int EqualScoringID; // Number of equal scoring result for this read
	uint EqualScoringCount; // Total number of equal scoring results
	//std::vector<LocationScore*> Scores;
	LocationScore * Scores;
	int nScores;
	int iScores;
	MappedRead * Paired; // Read pair
	uint Status;
	float Identity;
	int NM;
	char Strand;
	int QStart;
	int QEnd;
	int mappingQlty;
	float s;

	int length;

	char * RevSeq;
	char * Seq;

	char * qlty;

	uint nameLength;
	char * name;

	char * Buffer1;
	char * Buffer2;

	static volatile int maxSeqCount;

#ifdef INSTANCE_COUNTING
	static volatile int sInstanceCount;
#endif

	MappedRead(int const readid, int const qrymaxlen);
	~MappedRead();

	// Saves score to this Read and returns a pointer to the saved object
//	LocationScore * AddScore(LocationScore const & score);
	LocationScore * AddScore(CSTableEntry const & score);
	LocationScore * AddScore(float const score, uint const loc, bool const reverse);

	void clearScores(bool const keepTop = true);
	int numScores() const;

	LocationScore * TLS() const;
	bool hasCandidates() const;
//	void TopN(uint n);

	void AllocBuffers();
	void AllocScores(LocationScore * tmp, int const n);

	inline void SetFlag(ReadStatus const flag) {
		Status |= flag;
	}

	inline bool HasFlag(ReadStatus const flag) const {
		return (Status & flag) != 0;
	}

	char const * computeReverseSeq();

private:

	//void Unique();



	void DeleteReadSeq();

	static bool SortPred(LocationScore * lhs, LocationScore * rhs) {
		return lhs->Score.f > rhs->Score.f;
	}
	static bool UniquePred(LocationScore * lhs, LocationScore * rhs) {
		if (lhs == 0 || rhs == 0)
			return false;
		return (lhs->Location.m_Location == rhs->Location.m_Location)
				&& (lhs->Location.m_RefId == rhs->Location.m_RefId);
	}
	static bool IsZero(LocationScore * arg) {
		return (arg == 0);
	}
};

#endif
