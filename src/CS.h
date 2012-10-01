#ifndef __CS_H__
#define __CS_H__

#include "NGM.h"

#include "SequenceLocation.h"
#include "MappedRead.h"
#include "RefEntry.h"

#include <map>
#include <stdlib.h>
#include <string>

#include <cmath>

//#include "CSCache.h"

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
	int m_WrittenReads;
	int m_DiscardedReads;
	//	volatile int m_Candidates;
	int m_Mode;
	bool m_EnableBS;
	int m_Overflows;
	//float weightSum;
	const IRefProvider* m_RefProvider;
	RefEntry* m_entry;
	int c_SrchTableBitLen;
	int c_BitShift;
	int c_SrchTableLen;
	uint m_PrefixBaseSkip;
	bool m_Fallback;
	typedef void (*PrefixIterationFn)(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	typedef void (CS::*AddLocationFn)(const SequenceLocation& loc, const double freq);
	static void BuildPrefixTable(ulong prefix, uint pos, void* data);
	static void PrefixSearch(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);

	static void PrefixMutateSearchEx(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data, int mpos = 0);
	virtual int CollectResultsStd(MappedRead* read);
	int CollectResultsFallback(MappedRead* read);
	void FilterScore(LocationScore* score);
	void CheckFallback();
	virtual void RunBatch();
	void SendToBuffer(MappedRead* read);
	CSTableEntry* rTable; // standard
	int currentState;
	int* rList;
	int rListLength;
	float m_CsSensitivity;
	float currentThresh;
	float maxHitNumber;

	inline uint Hash(uint n) {
		//Multiplication Method (Corment)
		//		static float A = 0.5f * (sqrt(5) - 1);
		//		static uint m = floor(A * pow(2, 32));
		static uint m = 2654435761;
		return ((n) * m) >> c_BitShift;
	}

	inline uint GetBin(uint pos) {
		//static int shift = calc_binshift(Config.GetInt("corridor"));
		static int shift = calc_binshift(12);
		return pos >> shift;
	}

	inline uint ResolveBin(uint bin) {
//		static int shift = calc_binshift(Config.GetInt("corridor"));
		static int shift = calc_binshift(12);
		static uint offset = (shift > 0) ? 1 << (shift - 1) : 0;
		return (bin << shift) + offset;
	}

private:
	static const int estimateCount = 40000;
	void debugCS(MappedRead * read, int& n, float& mi_Threshhold);

	LocationScore * tmp;
	int tmpSize;

public:
	inline static int calc_binshift(int corridor) {
		corridor >>= 1;
		int l = 0;
		while ((corridor >>= 1) > 0)
			++l;
		return l;
	}
	//AddLocationFn AddLocation;
	virtual void AddLocationStd(const uint loc, const bool reverse, const double freq);
	void AddLocationFallback(const SequenceLocation& loc, const double freq);
	static uint prefixBasecount;
	static uint prefixBits;
	static ulong prefixMask;

	inline const int GetMode() {
		return m_Mode;
	}

	inline static const uint GetPrefixLength() {
		return prefixBasecount;
	}

	static void Init();
	static void Cleanup();
	static void PrefixMutateSearch(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void PrefixIteration(const char* sequence, uint length, PrefixIterationFn func, ulong mutateFrom, ulong mutateTo, void* data, uint prefixskip = 0,
			uint offset = 0);
	static void PrefixIteration(const char* sequence, uint length, PrefixIterationFn func, ulong mutateFrom, ulong mutateTo, void* data, uint prefixskip, uint offset,
			int prefixBaseCount);
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

//#define CallMemberFn(obj,pFn)  ((obj).*(pFn))

#endif
