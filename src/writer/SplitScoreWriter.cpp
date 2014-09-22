
/*
 * SplitScoreWriter.cpp
 *
 *  Created on: Aug 7, 2013
 *      Author: philipp_
 */

#include "SplitScoreWriter.h"

#include <string.h>

#include "Config.h"
#include "SequenceProvider.h"
#include "NGM.h"

//Format: http://samtools.sourceforge.net/SAM1.pdf
static int const report_offset = 1;

void SplitScoreWriter::DoWriteProlog() {

	char const * refName = 0;
	int refNameLength = 0;

	Print("#%u\n", NGM.Stats->TotalSeqs);
	Print("#");
	for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
		refName = SequenceProvider.GetRefName(i, refNameLength);
		Print("%d:%.*s\t", i/2, refNameLength, refName);
		if (NGM.DualStrand())
			++i;
		m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer, false);
	}
	Print("\n");
	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer, true);
}

void SplitScoreWriter::DoWriteRead(MappedRead const * const read, int const scoreID) {
	DoWriteReadGeneric(read, scoreID, "*", -1, 0, read->mappingQlty);
	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer);
}

void SplitScoreWriter::DoWriteReadGeneric(MappedRead const * const read, int const scoreID, char const * pRefName, int const pLoc, int const pDist, int const mappingQlty, int flags) {

//	Log.Message("Writing read %s", read->name);
	int n = 0;
	while (n < read->numScores() && read->Scores[n].Score.f > 0.0f) {
		n += 1;
	}

	Print("%s", read->name);
	int last = 0;
	for (int j = 0; j < n; ++j) {
		//fprintf(ofp2, "\t%u:%f", cur_read->Scores[j].Location.m_Location, cur_read->Scores[j].Score.f);

		int current = (int) read->Scores[j].Score.f;
		Print("\t%d:%u:%d:%d", read->Scores[j].Location.getrefId() / 2, read->Scores[j].Location.m_Location, read->Scores[j].Location.isReverse(), abs(last - current));
		last = current;
	}
	Print("\n");

//	Print("MD:Z:%s", read->Alignments[scoreID].pBuffer2);
//
//	Print("\n");

}

void SplitScoreWriter::DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2) {
	throw "Not implemented";
}

void SplitScoreWriter::DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, char const pRefName, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags =
		0) {
	//SRR002320.10000027.1    4       *       0       0       *       *       0       0       TTTATGTTGTTAATGTGTTGGGTGAGTGCGCCCCAT    IIIIIIIIIIIIIIIIIIIIIIIIII

	NGM.AddUnmappedRead(read, MFAIL_NOCAND);

}

void SplitScoreWriter::DoWriteUnmappedRead(MappedRead const * const read, int flags) {
	DoWriteUnmappedReadGeneric(read, -1, '*', -1, -1, 0, 0, flags | 0x04);
	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer);
}

void ScoreWriter::DoWriteEpilog() {

}


