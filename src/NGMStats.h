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
