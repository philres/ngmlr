#ifndef __FFORMATWRITER_H__
#define __FFORMATWRITER_H__

#include "GenericReadWriter.h"

class FFormatWriter : public GenericReadWriter
{
public:
	FFormatWriter(char const * const filename) : GenericReadWriter( filename )
	{
	}

protected:
	virtual void DoWriteProlog();
	virtual void DoWriteRead(MappedRead const * const read);
	virtual void DoWritePair(MappedRead const * const read1, MappedRead const * const read2);
	virtual void DoWriteUnmappedRead(
				MappedRead const * const read);
	virtual void DoWriteEpilog();
};

#endif
