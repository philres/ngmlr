/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
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

private:

public:
	void WriteProlog() {
		DoWriteProlog();
	}

	void WriteRead(MappedRead const * const read, bool mapped = true) {

		if (mapped) {
			std::map<SequenceLocation, bool> iTable;
			bool mappedOnce = false;
			for (int i = 0; i < read->Calculated; ++i) {


				float minResidues = Config.getMinResidues();

				if (minResidues <= 1.0f) {
					minResidues = read->length * minResidues;
				}

				mapped = mapped && (read->Alignments[i].Identity >= identity);
				mapped = mapped && ((float) (read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd) >= minResidues);

				Log.Debug(4, "READ_%d\tOUTPUT\tChecking alignment CRM_%d (read length: %d - %d - %d)\t%f >= %f\t%f >= %f", read->ReadId, i, read->length, read->Alignments[i].QStart, read->Alignments[i].QEnd, read->Alignments[i].Identity, identity, (float)(read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd), minResidues);

				if (mapped) {
					mappedOnce = true;
					if (iTable.find(read->Scores[i].Location) == iTable.end()) {
						iTable[read->Scores[i].Location] = true;
						Log.Debug(4, "READ_%d\tOUTPUT\tWriting alignment CRM_%d", read->ReadId, i);
						static const bool printAll = Config.getPrtintAll();
						if (!read->Alignments[i].skip || printAll) {
							DoWriteRead(read, i);
						}
					} else {
						Log.Debug(4, "READ_%d\tOUTPUT\tIgnoring duplicated alignment CRM_%d", read->ReadId, i);
					}
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
