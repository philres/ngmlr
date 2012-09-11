#ifndef __IREFPROVIDER_H__
#define __IREFPROVIDER_H__

#include "Types.h"

struct Location {
	uint m_Location;

	bool used() {
		return m_Location != 0;
	}
};

struct RefEntry {
	//static const int cIncrements[] = { 10, 20, 50, 100 };
	static int MaxRefsPerEntry;

	RefEntry(int locs = MaxRefsPerEntry) {
		reverse = false;
		weight = 0.0f;
		refCount = 0;
		refTotal = 0;
		nextEntry = 0;
	}
	~RefEntry() {
		if (nextEntry != 0)
			delete nextEntry;

		delete[] ref;
	}
	Location * ref;
	bool reverse;
	float weight;
	int refCount;
	int refTotal;
	RefEntry* nextEntry;
private:
	//---section Flattening
	RefEntry(void const * pEntries, int count);
	friend class PrefixTable;
	//---section Flattening

};

class IRefProvider {
public:

	virtual RefEntry const * GetRefEntry(ulong prefix,
			RefEntry * entry = 0) const = 0;

	virtual ~IRefProvider() {
	}
	;
};

#endif
