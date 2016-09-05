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

#ifndef __CIGARWRITER_H__
#define __CIGARWRITER_H__

#include "GenericReadWriter.h"

class SAMWriter: public GenericReadWriter {
public:

	SAMWriter(FileWriter * writer) :
			GenericReadWriter() {

//		if(Config.Exists(RG_ID)) {
//			RG = Config.GetString(RG_ID);
//		} else {
			RG = 0;
//		}
		m_Writer = writer;
	}

	virtual ~SAMWriter() {
		m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer, true);
	}

protected:
	virtual void DoWriteProlog();
	virtual void DoWriteRead(MappedRead const * const read, int const scoreID);
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2);
	virtual void DoWriteReadGeneric(MappedRead const * const read, int const scoreID, char const * pRefName, int const pLoc, int const pDist, int const mappingQlty, int flags =
			0);
	virtual void DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, char const pRefName, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags = 0x4);
	virtual void DoWriteEpilog();

private:

	FileWriter * m_Writer;

	char const * RG;
};

#endif
