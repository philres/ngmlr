#ifndef __CIGARWRITER_H__
#define __CIGARWRITER_H__

#include "GenericReadWriter.h"

class SAMWriter: public GenericReadWriter {
	public:
		SAMWriter(char const * const filename) :
			GenericReadWriter(filename) {

		}

	protected:
		virtual void DoWriteProlog();
		virtual void DoWriteRead(MappedRead const * const read);
		virtual void DoWritePair(MappedRead const * const read1, MappedRead const * const read2);
		virtual void DoWriteReadGeneric(MappedRead const * const read, char const * pRefName, int const pLoc, int const pDist, int const mappingQlty, int flags = 0);
		virtual void DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, char const pRefName, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags);
		virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags = 0x4);
		virtual void DoWriteEpilog();
};

#endif
