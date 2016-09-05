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

#ifndef __LOCATIONSCORE_H__
#define __LOCATIONSCORE_H__

#include "SequenceLocation.h"
#include "PlatformSpecifics.h"

struct MappedRead;

struct CSTableEntry {
	uint m_Location;
	uint state;
	float fScore;
	float rScore;
};

union UScore {
	int i;
	float f;
};

struct LocationScore {
	UScore Score;
	SequenceLocation Location;

#ifdef INSTANCE_COUNTING
	static volatile int sInstanceCount;
#endif

	LocationScore() {
#ifdef INSTANCE_COUNTING
		AtomicInc(&sInstanceCount);
#endif
	}

	~LocationScore() {
#ifdef INSTANCE_COUNTING
		AtomicDec(&sInstanceCount);
#endif
	}

	LocationScore(SequenceLocation location, int score) {
		Score.i = score;
		Location = location;

#ifdef INSTANCE_COUNTING
		AtomicInc(&sInstanceCount);
#endif
	}

};

#endif
