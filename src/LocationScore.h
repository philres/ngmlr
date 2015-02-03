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
