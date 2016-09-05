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

	long long corridorLen;

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

		corridorLen = 0;

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
