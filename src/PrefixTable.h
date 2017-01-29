/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef COMPACTPREFIXTABLE_H_
#define COMPACTPREFIXTABLE_H_

#include "IRefProvider.h"
#include "Types.h"

#include <map>
#include <vector>

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

//Represents a single hashtable unit
//For genomes < 4GB, we always only use one unit with no genomic offset
//In order to support larger genomes (genomic kmer position can be > INT_MAX),
//we create another hashtable unit, assign it an offset value, and store all the large
//kmer positions there.
struct TableUnit
{
	 TableUnit() : cRefTableLen(0), RefTable(0), RefTableIndex(0) { Offset = 0; }
	~TableUnit()
	{
		delete[] RefTable;
		RefTable = 0;

		delete[] RefTableIndex;
		RefTableIndex = 0;
	}

	uint      cRefTableLen;
	Location* RefTable;
	Index*    RefTableIndex;

	uloc      Offset;
};

#pragma pack(pop)

class CompactPrefixTable: public IRefProvider {
public:
	CompactPrefixTable(bool const dualStrand = true, bool const skip = true);
	virtual ~CompactPrefixTable();

	const RefEntry* GetRefEntry(ulong prefix, RefEntry* entry) const;
	uint GetRefEntryChainLength() const;

	static int maxPrefixFreq;

private:
	//Biggest location value supported by a single table (precision of used location data type)
	//=> Determines: How many table units do we need?
	//   Table units created:
	//      Reference genome size divided by c_tableLocMax, offsets increasing every table by c_tableLocMax
	static uloc c_tableLocMax; //4294967296; //UINT_MAX

	//Table units array
	TableUnit* m_Units;
	uint m_UnitCount;

	//Static members used to direct generation of table units
	static TableUnit* CurrentUnit;


	/********************************************************************/
	/*************Used only for building reference table*****************/
	/********************************************************************/
	//Used to control which kmers should be counted for index building, only locations
	//that will be in the unit should also be in the index
	static uloc kmerCountMinLocation;
	static uloc kmerCountMaxLocation;

	int m_CurGenSeq;
	static ulong lastPrefix;
	static loc lastBin;
	static loc lastPos;
	static uint skipCount;
	static uint skipBuild;

	uint m_RefSkip;
	uint m_PrefixLength;
	/********************************************************************/

	//Config values
	bool DualStrand;
	bool skipRep;

	void Generate();
	
	void CreateTable(const uint length);
	int* CountKmerFreq(const uint length);
	uint createRefTableIndex(const uint length);
	static void CountKmer(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void CountKmerwoSkip(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void BuildPrefixTable(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo, void* data);
	static void BuildPrefixTablewoSkip(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo, void* data);
	void SaveToRefTable(ulong prefix, Location loc);
	void Clear();
	void saveToFile(const char* fileName, const uint refIndexSize);
	bool readFromFile(const char* fileName);

	void test();
};

#endif /* COMPACTPREFIXTABLE_H_ */
