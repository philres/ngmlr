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
	void DecodeRefSequence(char* const buffer, int n, uint offset, uint len);
	// Gets the length of read/reference string n
	virtual uint GetRefLen(int n) const;
	virtual uint GetConcatRefLen() const;
	virtual int GetRefCount() const;
	const virtual char* GetRefName(int n, int& len) const;
	virtual uint GetRefStart(int n) const;

	virtual void PagingUpdate();
	static const int maxRefNameLength = 100;

	SequenceLocation convert(MappedRead * read, uint m_Location);

private:

	//static const size_t MAX_REF_NAME_LENGTH;

	_SequenceProvider();
	~_SequenceProvider();
//	NGMMutex f_mutex;
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

	//int refIdxMap;
	//int refMap;
	//int refcplMap;
	//int refBaseMap;

	int refCount;

	bool m_EnableBS;

	void writeEncRefToFile(char const * fileName, uint const refCount, uint const encRefSize);
	int readEncRefFromFile(char const * fileName);

//	void PrepareInputFile(std::string);
//	void readFromPosition(FILE* file, char* str, ulong offset,
//			const int len) const;
};

#define SequenceProvider _SequenceProvider::Instance()

#endif
