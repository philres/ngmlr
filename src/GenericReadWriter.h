/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __GENERICREADWRITER_H__
#define __GENERICREADWRITER_H__

#include <map>
#include <stdio.h>
#include <stdarg.h>

#include "ILog.h"
#include "IConfig.h"
#include "NGM.h"
#include "MappedRead.h"
#include "FileWriter.h"

#undef module_name
#define module_name "OUTPUT"

class GenericReadWriter {

protected:

	virtual void DoWriteProlog() = 0;
	virtual void DoWriteRead(MappedRead const * const read,
			int const scoreID) = 0;
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1,
			MappedRead const * const read2, int const scoreId2) = 0;
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags =
			0x4) = 0;
	virtual void DoWriteEpilog() = 0;

	float identity;

	bool writeUnmapped;

	static int const BUFFER_SIZE =  3200000;
	static int const BUFFER_LIMIT = 1600000;

	char * writeBuffer;
	int bufferPosition;

	int Print(const char *format, ...) {
		int done;
		va_list arg;

		va_start(arg, format);
		done = vsprintf(writeBuffer + bufferPosition, format, arg);
		bufferPosition += done;
		if(bufferPosition >= BUFFER_SIZE) {
			Log.Error("Size of write buffer exceeded");
		}
		va_end(arg);
		return done;
	}

public:

	GenericReadWriter() {
		writeBuffer = new char[BUFFER_SIZE];
		bufferPosition = 0;

		identity = Config.getMinIdentity();
		writeUnmapped = Config.getWriteUnampped();
	}
	virtual ~GenericReadWriter() {
		if (writeBuffer != 0) {
			delete[] writeBuffer;
			writeBuffer = 0;
		}
	}

	void WriteProlog() {
		DoWriteProlog();
	}

	void WriteRead(MappedRead const * const read, bool mapped = true) {

		bool mappedOnce = false;
		for (int i = 0; i < read->Calculated && mapped; ++i) {
			mappedOnce = mappedOnce || !read->Alignments[i].skip;
			if(!read->Alignments[i].skip) {
				DoWriteRead(read, i);
			}
		}
		if (mappedOnce) {
			Log.Debug(4, "READ_%d\tOUTPUT\tRead was mapped", read->ReadId);
			NGM.AddMappedRead(read->ReadId);
		} else {
			if(read->HasFlag(NGMNames::Empty)) {
				Log.Debug(4, "READ_%d\tOUTPUT\tRead empty (discard read)", read->ReadId);
			} else {
				Log.Debug(4, "READ_%d\tOUTPUT\tRead unmapped", read->ReadId);
				DoWriteUnmappedRead(read);
			}
		}
	}

	void WriteEpilog() {
		DoWriteEpilog();
	}
}
;

#endif
