/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __SEQUENCEPROVIDER_H__
#define __SEQUENCEPROVIDER_H__

#include <string>
#include <map>

#include "Types.h"
#include "NGMThreads.h"
#include "SequenceLocation.h"
#include "MappedRead.h"

class _SequenceProvider {
public:
	static _SequenceProvider & Instance();
	static void Cleanup();

	virtual void Init(bool dualstrand = true);

	struct Chromosome {
		uloc start;
		uloc end;
	};

	bool DecodeRefSequenceExact(char * const buffer, uloc offset,
			uloc bufferLength, int corridor);
	bool DecodeRefSequence(char* const buffer, int n, uloc offset, uloc len);
	// Gets the length of read/reference string n
	virtual uloc GetRefLen(int n) const;
	virtual uloc GetConcatRefLen() const;
	virtual int GetRefCount() const;
	const virtual char* GetRefName(int n, int& len) const;
	virtual uloc GetRefStart(int n) const;

	virtual void PagingUpdate();
	static const int maxRefNameLength = 100;

	bool convert(SequenceLocation & m_Location);

	Chromosome getChrBorders(loc start, loc stop);
	Chromosome getChrStart(uloc const position);

private:

	_SequenceProvider();
	~_SequenceProvider();

	bool CheckQryNr(int n) const;
	bool CheckRefNr(int n) const;
	bool DualStrand;
	static _SequenceProvider* pInstance;

	struct RefIdx {
		uint SeqId;
		uint Flags;
		uloc SeqStart;
		uint SeqLen;
		uint NameLen;
		char name[maxRefNameLength];
	};

	RefIdx * binRefIdx;
	char* binRef;
	uloc binRefIndex;

	// Files
	std::string refFileName;
	std::string refBaseFileName;

	uloc binRefSize;

	int refCount;

	uloc * refStartPos;

	static const int minRefSeqLen = 10;

	void writeEncRefToFile(char const * fileName, uint const refCount,
			uloc const encRefSize);
	int readEncRefFromFile(char const * fileName, const uloc maxLen);
	uloc decode(uloc startPosition, uloc endPosition, char* const sequence);
};

#define SequenceProvider _SequenceProvider::Instance()

#endif
