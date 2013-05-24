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
	MappedRead * Read;

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

	LocationScore(SequenceLocation location, int score, MappedRead * read) {
		Score.i = score;
		Location = location;
		Read = read;

#ifdef INSTANCE_COUNTING
		AtomicInc(&sInstanceCount);
#endif
	}


//	LocationScore(uint const loc, bool const reverse, float const score, MappedRead * read) {
//		Score.f = score;
//		Location.m_Location = loc;
//		Location.m_RefId = reverse;
//		Read = read;
//#ifdef INSTANCE_COUNTING
//		AtomicInc(&sInstanceCount);
//#endif
//	}
};

#endif
