#ifndef __IREFPROVIDER_H__
#define __IREFPROVIDER_H__

#include "Types.h"

struct Location {
	uint m_Location;

	bool used() {
		return m_Location != 0; //TODO_GENOMESIZE: May be null (first loc in second hashtable/)
	}
};

struct RefEntry {
	//static const int cIncrements[] = { 10, 20, 50, 100 };
	//static int MaxRefsPerEntry; //UNUSED

	RefEntry(int locs = 0) {
		reverse = false;
		weight = 0.0f;
		refCount = 0;
		refTotal = 0;
	}
	~RefEntry() {
		//delete[] ref;
	}

	inline uloc getRealLocation(const Location& loc) const { return loc.m_Location + offset; }

	Location * ref;
	uloc offset; //Offset for all Locations
	bool reverse;
	float weight;
	int refCount;
	int refTotal;
private:
	RefEntry(void const * pEntries, int count);
	friend class PrefixTable;

};

class IRefProvider {
public:

	virtual RefEntry const * GetRefEntry(ulong prefix,
			RefEntry * entry = 0) const = 0;

	virtual uint GetRefEntryChainLength() const = 0;

	virtual ~IRefProvider() {
	}
	;
};

#endif
