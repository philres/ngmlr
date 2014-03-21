/*
 * ScoreWriter.h
 *
 *  Created on: Aug 7, 2013
 *      Author: philipp_
 */

#ifndef SCOREWRITER_H_
#define SCOREWRITER_H_

#include "GenericReadWriter.h"
#include "FileWriter.h"

class ScoreWriter: public GenericReadWriter {
public:
//		SAMWriter(char const * const filename) :
//			GenericReadWriter(filename) {
	ScoreWriter(FileWriter * writer) :
	GenericReadWriter() {

		m_Writer = writer;
	}

	virtual ~ScoreWriter() {
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
};

#endif /* SCOREWRITER_H_ */
