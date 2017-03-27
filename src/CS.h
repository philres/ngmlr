/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __CS_H__
#define __CS_H__

//#include "NGM.h"

#include "SequenceLocation.h"
#include "MappedRead.h"
#include "ScoreBuffer.h"

#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

#include <cmath>

/*
 Candidate Search:

 Initialisierung:
 Iteriere ueber jeden Praefix der definierten Referenzen (config:ref_mode mit -1=alle, 0..n eine einzelne referenz)
 Fuege in RefTable fuer diesen Praefix seine Position und die Id des Referenzstrings ein

 Suche:
 Iteriere ueber jeden Pruefix des Reads
 Lies fuer diesen Praefix die Positionsliste aus RefTable
 Fuer jede Position dieses Praefixes
 Offset um die Position des Praefixes (=normalisierung auf den Read-Anfang)
 fast_mode:
 Pruefe rTable-Capactiy, wenn overflow -> Neustart der Suche im safe_mode
 Schreibe in rTable mit linearem Kollisionshandling
 save_mode:
 Schreibe in iTable

 Sammle Ergebnisse ein:
 Iteriere ueber jedes Element in fast_mode:rTable / save_mode:iTable
 Wenn count > threshhold, schicke es zur weiterverarbeitung
 */

class CS: public NGMTask {
protected:
	static volatile int s_ThreadCount;
	const int m_CSThreadID;
	const int m_BatchSize;
	std::vector<MappedRead*> m_CurrentBatch;
	int m_CurrentSeq;
	int m_CurrentReadLength;
	int m_ProcessedReads;
	int m_DiscardedReads;
	//	volatile int m_Candidates;

	int m_Overflows;
	//float weightSum;
	const IRefProvider* m_RefProvider;
	RefEntry* m_entry;
	uint m_entryCount;
	int c_SrchTableBitLen;
	int c_BitShift;
	int c_SrchTableLen;
	uint m_PrefixBaseSkip;
	bool m_Fallback;
	typedef void (*PrefixIterationFn)(ulong prefix, uloc pos, ulong mutateFrom,
			ulong mutateTo, void* data);
	typedef void (CS::*AddLocationFn)(const SequenceLocation& loc,
			const double freq);
	static void BuildPrefixTable(ulong prefix, uloc pos, void* data);
	static void PrefixSearch(ulong prefix, uloc pos, ulong mutateFrom,
			ulong mutateTo, void* data);

	virtual int CollectResultsStd(MappedRead* read);
	int CollectResultsFallback(MappedRead* read);
	void FilterScore(LocationScore* score);
	void CheckFallback();
	virtual int RunBatch(ScoreBuffer * sw, AlignmentBuffer * out);
	int RunRead(MappedRead * currentRead, PrefixIterationFn pFunc,
			ScoreBuffer * sw, AlignmentBuffer * out);
	void SendToBuffer(MappedRead* read, ScoreBuffer * sw,
			AlignmentBuffer * out);

	void AllocRefEntryChain();

	CSTableEntry* rTable; // standard
	int currentState;
	int* rList;
	int rListLength;
	float m_CsSensitivity;
	float currentThresh;
	float maxHitNumber;

	uint hpoc;

	inline uint Hash(uloc n) {
		//Multiplication Method (Corment)
		//static float A = 0.5f * (sqrt(5) - 1);
		//static uloc m = floor(A * pow(2, 64));
		//static uint m = 2654435761;
		static uloc m = 11400714819323199488u;

		return uint((n * m) >> c_BitShift);
	}

	inline void SetSearchTableBitLen(int bitLen) {
		if (bitLen >= 24) {
			Log.Error("SearchTable exceeded length.");
		}
		c_SrchTableBitLen = bitLen;
		c_BitShift = 64 - c_SrchTableBitLen;
		c_SrchTableLen = (int) pow(2, c_SrchTableBitLen);
	}

private:

	static const int estimateCount = 40000;
	void debugCS(MappedRead * read, int& n, float& mi_Threshhold);
	Align computeAlignment(MappedRead* read, int const scoreId,
			int const corridor);

	float scoreReadPart(char const * const readSeq, int const qryLen, int const bin,
			long const refBin);

	LocationScore * tmp;
	int tmpSize;
	IAlignment * oclAligner;
	AlignmentBuffer * alignmentBuffer;

public:

	//AddLocationFn AddLocation;
	virtual void AddLocationStd(const uloc loc, const bool reverse,
			const double freq);
	void AddLocationFallback(const SequenceLocation& loc, const double freq);
	static uint prefixBasecount;
	static uint prefixBits;
	static ulong prefixMask;

	inline static const uint GetPrefixLength() {
		return prefixBasecount;
	}

	static void Init();
	static void Cleanup();
	static void PrefixIteration(const char* sequence, uloc length,
			PrefixIterationFn func, ulong mutateFrom, ulong mutateTo,
			void* data, uint prefixskip = 0, uloc offset = 0);
	static void PrefixIteration(const char* sequence, uloc length,
			PrefixIterationFn func, ulong mutateFrom, ulong mutateTo,
			void* data, uint prefixskip, uloc offset, int prefixBaseCount);
	CS(bool useBuffer = true);
	~CS();
	void DoRun();

	inline int GetStage() const {
		return 0;
	}

	inline const char* GetName() const {
		return "CS";
	}
};

//int const static csBitShift = 2;

//inline static int calc_binshift() {
//	int corridor = 4;
//	int l = 0;
//	while ((corridor >>= 1) > 0)
//		++l;
//	return l;
//}

inline uloc GetBin(uloc pos) {
	static int const csBitShift = Config.getBinSize();
//	static int shift = calc_binshift();
	return pos >> csBitShift;
}

inline uloc ResolveBin(uloc bin) {
	static int const csBitShift = Config.getBinSize();
//	static int shift = calc_binshift();
	static uint offset = (csBitShift > 0) ? 1 << (csBitShift - 1) : 0;
	return (bin << csBitShift) + offset;
}

#endif
