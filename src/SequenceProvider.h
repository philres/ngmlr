#ifndef __SEQUENCEPROVIDER_H__
#define __SEQUENCEPROVIDER_H__

#include <string>
#include <map>

#include "Types.h"
#include "NGMThreads.h"

#include "SequenceLocation.h"
#include "MappedRead.h"

/*
 *  SequenceProvider: Schnittstelle zu den Sequenz-Daten
 */
class _SequenceProvider {
public:
	static _SequenceProvider & Instance();
	static void Cleanup();

	virtual void Init();
	bool DecodeRefSequence(char* const buffer, int n, uint offset, uint len);
	// Gets the length of read/reference string n
	virtual uint GetRefLen(int n) const;
	virtual uint GetConcatRefLen() const;
	virtual int GetRefCount() const;
	const virtual char* GetRefName(int n, int& len) const;
	virtual uint GetRefStart(int n) const;

	virtual void PagingUpdate();
	static const int maxRefNameLength = 100;

	bool convert(SequenceLocation & m_Location);

private:

	_SequenceProvider();
	~_SequenceProvider();

	bool CheckQryNr(int n) const;
	bool CheckRefNr(int n) const;
	static _SequenceProvider* pInstance;

	struct RefIdx {
		uint SeqId;
		uint Flags;
		ulong SeqStart;
		uint SeqLen;
		uint NameLen;
		char name[maxRefNameLength];
	};

	RefIdx * binRefIdx;
	char* binRef;
	int binRefIndex;

	// Files
	std::string refFileName;
	std::string refBaseFileName;

	long refFileLen;
	long refCplFileLen;

	long refBaseFileLen;

	uint binRefSize;

	int refCount;

	bool m_EnableBS;

	uint * refStartPos;

	void writeEncRefToFile(char const * fileName, uint const refCount, uint const encRefSize);
	int readEncRefFromFile(char const * fileName);

};

#define SequenceProvider _SequenceProvider::Instance()

#endif
