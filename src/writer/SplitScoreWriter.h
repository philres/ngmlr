/*
 * ScoreWriter.h
 *
 *  Created on: Aug 7, 2013
 *      Author: philipp_
 */

#ifndef SPLITSCOREWRITER_H_
#define SPLITSCOREWRITER_H_

#include "GenericReadWriter.h"
#include "FileWriter.h"
#include "zlib.h"

class SplitScoreWriter: public GenericReadWriter {
public:

	gzFile m_Output_1;
	gzFile m_Output_2;
	gzFile m_Output_3;
	gzFile m_Output_4;
	gzFile m_Output_5;


	SplitScoreWriter(char const * const filename) :
			GenericReadWriter() {

	}

	virtual ~SplitScoreWriter() {

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


};

#endif /* SPLITSCOREWRITER_H_ */
