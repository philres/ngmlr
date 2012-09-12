#ifndef ___BAMWRITER_H___
#define ___BAMWRITER_H___

#include "GenericReadWriter.h"

#ifdef _BAM

#include "api/BamAlignment.h"

class BAMWriter: public GenericReadWriter {
public:
	BAMWriter(char const * const filename) :
			GenericReadWriter(filename), file(filename) {
	}

protected:
	virtual void DoWriteProlog();
	virtual void DoWriteRead(MappedRead const * const read);
	virtual void DoWritePair(MappedRead const * const read1,
			MappedRead const * const read2);
	virtual void DoWriteReadGeneric(MappedRead const * const read,
			char const pRefName, int const pLoc, int const pDist,
			int const mappingQlty, int flags = 0);
	virtual void DoWriteUnmappedReadGeneric(MappedRead const * const read,
			int const refId, char const pRefName, int const loc, int const pLoc,
			int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedRead(MappedRead const * const read);
	virtual void DoWriteEpilog();

private:
	char const * const file;
	void translate_flag(BamTools::BamAlignment &al, int flags);

};
#endif

#endif //___BAMWRITER_H___
