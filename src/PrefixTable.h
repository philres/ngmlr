/*
 * CompactPrefixTable.h
 *
 *  Created on: May 9, 2012
 *      Author: philipp_
 */

#ifndef COMPACTPREFIXTABLE_H_
#define COMPACTPREFIXTABLE_H_

#include "IRefProvider.h"
#include "Types.h"

#pragma pack(push)
#pragma pack(1)

struct Index {

	uint m_TabIndex;
	char m_RevCompIndex;

	Index() {
		m_TabIndex = 0;
		m_RevCompIndex = 0;
	}

	bool used() {
		//return m_RevCompIndex;
		return m_RevCompIndex != 0;
		//return m_RevCompIndex & 0x80000000;
	}

//	void setUsed() {
//		//m_RevCompIndex = m_RevCompIndex | 0x80000000;
//		m_RevCompIndex = true;
//	}

//	void setRevCompIndex(uint prefix) {
//		if(used()) {
//			m_RevCompIndex = prefix;
//			setUsed();
//		} else {
//			m_RevCompIndex = prefix;
//		}
//	}
//
//	uint getRevCompIndex() {
//		return m_RevCompIndex & 0x7FFFFFFF;
//	}

};

#pragma pack(pop)

class CompactPrefixTable: public IRefProvider {
public:
	CompactPrefixTable(bool const dualStrand = true, bool const skip = true);
	virtual ~CompactPrefixTable();

	const RefEntry* GetRefEntry(ulong prefix, RefEntry* entry) const;

	static int maxPrefixFreq;

private:

	uint cRefTableLen;
	Location* RefTable;
	Index* RefTableIndex;
	int m_CurGenSeq;
	int m_RECount;
	int m_RRCount;
	int m_BCalls;
	ulong m_TotalLocs;
	static ulong lastPrefix;
	static int lastBin;
	static int lastPos;
	static uint skipCount;
	static uint skipBuild;
	uint m_RefSkip;
	uint m_PrefixLength;
	bool DualStrand;
	bool skipRep;
	void Generate();
	void CreateTable(const uint length);
	int* CountKmerFreq(const uint length);
	uint createRefTableIndex(const uint length);
	static void CountKmer(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void CountKmerwoSkip(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void BuildPrefixTable(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void BuildPrefixTablewoSkip(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data);
	void SaveToRefTable(ulong prefix, Location loc);
	void Clear();
	void saveToFile(const char* fileName, const uint refIndexSize, const uint refTableSize);
	uint readFromFile(const char* fileName);

	ulong CalcSize() const;
	static bool CheckHeader(const char* const buffer);

};

#endif /* COMPACTPREFIXTABLE_H_ */
