#ifndef __NGMSTATS_H__
#define __NGMSTATS_H__

#include "Types.h"

struct NGMStats
{
	int TotalSeqs;

	float csTime;
	float scoreTime;
	float alignTime;
	int readsPerSecond;
	int csLength;
	int csOverflows;
	int avgnCRMS;

	float validPairs;
	float insertSize;

	int readsInProcess;

	int readObjectsInProcess;

	uint invalidAligmentCount;
	uint alignmentCount;

public:
	NGMStats() {
		TotalSeqs = 0;
		csTime = 0.0f;
		scoreTime = 0.0f;
		alignTime = 0.0f;
		readsPerSecond = 0;
		csLength = 0;
		csOverflows = 0;
		avgnCRMS = 0;

		validPairs = 0.0f;
		insertSize = 0.0f;

		readsInProcess = 0;
		readObjectsInProcess = 0;

		invalidAligmentCount = 0;
		alignmentCount = 0;
	}

	~NGMStats() {

	}

};

#endif
