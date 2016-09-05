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

#ifndef ___BAMWRITER_H___
#define ___BAMWRITER_H___



#include "GenericReadWriter.h"
#include "FileWriterBam.h"
#include <string.h>

class BAMWriter: public GenericReadWriter {
public:
//	BAMWriter(char const * const filename) :
//			GenericReadWriter(filename), file(filename) {
	BAMWriter(FileWriterBam * pWriter, char const * const pFile) :
			GenericReadWriter(), writer(pWriter), file(pFile) {
		NGMInitMutex(&m_OutputMutex);
		bufferIndex = 0;

//		if (Config.Exists(RG_ID)) {
//			RG = std::string(Config.GetString(RG_ID));
//		}
	}

	~BAMWriter() {
		if (bufferIndex > 0) {
			writer->SaveAlignment(buffer, bufferIndex);
			bufferIndex = 0;
		}
	}

protected:
	virtual void DoWriteProlog();
	virtual void DoWriteRead(MappedRead const * const read, int const scoreId);
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2);
	virtual void DoWriteReadGeneric(MappedRead const * const read, int const scoreId, int const pRef, int const pLoc, int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, int const pRef, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags = 0x4);
	virtual void DoWriteEpilog();

private:
	void translate_flag(BamTools::BamAlignment * al, int flags);
	void addAdditionalInfo(const MappedRead* const read, BamTools::BamAlignment* al);

	FileWriterBam * writer;

	char const * const file;

	NGMMutex m_OutputMutex;

	BamTools::BamAlignment * buffer[10000];

	int bufferIndex;

	std::string RG;

};


#endif //___BAMWRITER_H___
