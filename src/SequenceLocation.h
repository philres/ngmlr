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

#ifndef __SEQUENCELOCATION_H__
#define __SEQUENCELOCATION_H__

#include <stdlib.h>
#include <assert.h>

#include "Types.h"
#include "IRefProvider.h"

struct SequenceLocation {
	uloc m_Location;

	bool used() const {
		return m_Location != 0;
	}

	int getrefId() const {
		return m_RefId & 0x7FFFFFFF;
	}

	void setRefId(int const id) {
		//assert(id >= 0);
		m_RefId = id | (m_RefId & 0x80000000);
	}

	bool isReverse() const {
		//return m_Reverse;
		return m_RefId & 0x80000000;
	}

	void setReverse(bool const reverse) {
		if (reverse) {
			m_RefId |= 0x80000000;
		} else {
			m_RefId = getrefId();
		}
	}

	bool operator<(SequenceLocation const & rhs) const {
		if (m_Location < rhs.m_Location)
			return true;
		else if (m_Location == rhs.m_Location)
			//m_RefId instead of getRefId. Strand is important. Equal positions on the same reference but not on the same strand are not equal.
			//To objects are deemed equal if !(a < b) && !(b < a): important for map in ReadProvider
			return (m_RefId < rhs.m_RefId);

		return false;
	}

	SequenceLocation() {
		m_Location = 0;
		m_RefId = 0;
	}

	SequenceLocation(uloc const loc, short const refid, bool const reverse) {
		m_Location = loc;
		setRefId(refid);
		setReverse(reverse);
	}

	SequenceLocation(Location const & other, uloc offset) {
		m_Location = other.m_Location + offset;
		setRefId(0);
	}

private:
	//bool m_Reverse;
	//unsigned short m_RefId;
	unsigned int m_RefId;
};
//#pragma pack(pop)

#endif
