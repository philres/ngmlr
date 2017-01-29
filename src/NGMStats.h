/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __NGMSTATS_H__
#define __NGMSTATS_H__

#include "Types.h"

struct NGMStats
{
	float csTime;
	float scoreTime;
	float alignTime;
	int csLength;
	int csOverflows;
	int avgnCRMS;

	float avgAlignPerc;

	long long corridorLen;

	long long readLengthSum;

	int readsInProcess;

	uint invalidAligmentCount;
	uint alignmentCount;

public:
	NGMStats() {
		csTime = 0.0f;
		scoreTime = 0.0f;
		alignTime = 0.0f;
		csLength = 0;
		csOverflows = 0;
		avgnCRMS = 0;
		avgAlignPerc = 0.0f;
		corridorLen = 0;
		readLengthSum = 0ll;

		readsInProcess = 0;

		invalidAligmentCount = 0;
		alignmentCount = 0;
	}

	~NGMStats() {

	}

};

#endif
