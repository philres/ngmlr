#ifndef __SEQUENCELOCATION_H__
#define __SEQUENCELOCATION_H__

#include "Types.h"
#include "IRefProvider.h"
#include <stdlib.h>
#include <assert.h>


struct SequenceLocation {
	uint m_Location;

	bool used() const {
		return m_Location != 0;
	}

	int getrefId() const {
		return m_RefId & 0x7FFFFFFF;
	}

	void setRefId(int const id) {
		assert(id > 0);
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
			return (getrefId() < rhs.getrefId());

		return false;
	}

	SequenceLocation() {
		m_Location = 0;
		m_RefId = 0;
	}

	SequenceLocation(uint const loc, short const refid, bool const reverse) {
		m_Location = loc;
		setRefId(refid);
		setReverse(reverse);
	}

	SequenceLocation(Location const & other) {
		m_Location = other.m_Location;
		setRefId(0);
	}

private:
	//bool m_Reverse;
	//unsigned short m_RefId;
	unsigned int m_RefId;
};
//#pragma pack(pop)

#endif
