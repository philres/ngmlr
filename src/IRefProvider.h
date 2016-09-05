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
		ref = 0;
		offset = 0;
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
