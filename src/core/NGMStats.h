#ifndef __NGMSTATS_H__
#define __NGMSTATS_H__

#include "Types.h"

struct CSThreadStats
{
	uint Status;
	uint CurrentRead;
	uint MappedReads;
	uint UnmappedReads;
	uint FinishedBatches;
	uint TotalCandidates;

	float TotalTime;
};

struct NGMStats
{
	static const int cMaxThreadStats = 32;

	int Cookie;

	int CurrentPartition;
	int TotalPartitions;

	int CSStage;
	int TotalSeqs;
	int TotalBases;
	int CurrentSeq;
	int CurrentSeqBases;
	//int CurrentBase;
	//int CurrentTotal;
	int UniquePrefixCount;
	int RefEntryCount;

	int CachePageWriteCur;
	int CachePageWriteMax;

	int Buffer1;
	int Buffer2;

	int CSThreads;
	int SWThreads;

	float csTime;
	float scoreTime;
	float alignTime;
	int readsPerSecond;
	int csLength;
	int csOverflows;
	int avgnCRMS;

	float validPairs;
	float insertSize;

	CSThreadStats CS[cMaxThreadStats];

	static NGMStats * InitStats(char const * const AppName);
private:
	NGMStats();
	~NGMStats();
};

#endif
