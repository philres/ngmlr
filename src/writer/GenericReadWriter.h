#ifndef __GENERICREADWRITER_H__
#define __GENERICREADWRITER_H__

#include <map>

#include "ILog.h"
#include "Config.h"
#include "NGM.h"

#include "MappedRead.h"
#include "FileWriter.h"

class GenericReadWriter {

public:

	GenericReadWriter() {
		writeBuffer = new char[BUFFER_SIZE];
		bufferPosition = 0;

		identity = 0;
		if (Config.Exists("identity_tresh")) {
			identity = Config.GetFloat("identity_tresh", 0, 100);
			identity = identity / 100.0f;
		}
		writeUnmapped = !Config.GetInt("no_unal");
	}
	virtual ~GenericReadWriter() {

	}
protected:

	virtual void DoWriteProlog() = 0;
	virtual void DoWriteRead(MappedRead const * const read, int const scoreID) = 0;
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2) = 0;
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags = 0x4) = 0;
	virtual void DoWriteEpilog() = 0;

	float identity;

	bool writeUnmapped;

	static int const BUFFER_SIZE = 17000000;
	static int const BUFFER_LIMIT = 16000000;



	char * writeBuffer;
	int bufferPosition;

	int Print(const char *format, ...) {
		int done;
//		NGMLock(&m_OutputMutex);
		va_list arg;

		va_start(arg, format);
		done = vsprintf(writeBuffer + bufferPosition, format, arg);
		bufferPosition += done;
		va_end(arg);
//		NGMUnlock(&m_OutputMutex);
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

				static float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
				static float minResidues = Config.GetFloat("min_residues", 0, 1000);

				if (minResidues <= 1.0f) {
					minResidues = read->length * minResidues;
				}

				mapped = mapped && (read->Alignments[i].Identity >= minIdentity);
				mapped = mapped && ((read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd) >= minResidues);

				if (mapped) {
					mappedOnce = true;
					if(iTable.find(read->Scores[i].Location) == iTable.end()) {
						iTable[read->Scores[i].Location] = true;
						DoWriteRead(read, i);
					} else {
						Log.Message("Ignoring duplicated alignment %d for read %s.", i, read->name);
					}
				}
			}
			if (mappedOnce) {
				NGM.AddMappedRead(read->ReadId);
			} else {
				DoWriteUnmappedRead(read);
			}
		} else {
			DoWriteUnmappedRead(read);
		}
	}

	void WritePair(MappedRead * const read1, int const scoreId1, MappedRead * const read2, int const scoreId2) {

		//TODO: fix paired end!!! MULTI MAP!

		static float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
		static float minResidues = Config.GetFloat("min_residues", 0, 1000);

		if (minResidues < 1.0f) {
			minResidues = std::min(read1->length, read2->length) * minResidues;
		}

		bool mapped1 = read1->hasCandidates() && (read1->Alignments[scoreId1].Identity >= minIdentity)
		&& ((read1->length - read1->Alignments[scoreId1].QStart - read1->Alignments[scoreId1].QEnd) >= minResidues);
		bool mapped2 = read2->hasCandidates() && (read2->Alignments[scoreId2].Identity >= minIdentity)
		&& ((read2->length - read2->Alignments[scoreId2].QStart - read2->Alignments[scoreId2].QEnd) >= minResidues);

		if (!mapped1) {
			read1->clearScores();
			//NGM.AddUnmappedRead(read1, MFAIL_IDENT);
		} else {
			NGM.AddMappedRead(read1->ReadId);
		}
		if (!mapped2) {
			read2->clearScores();
			//NGM.AddUnmappedRead(read2, MFAIL_IDENT);
		} else {
			NGM.AddMappedRead(read2->ReadId);

		}

		//Log.Message("Output paired 1: hC %d, R: %d %d, I: %f %f", read1->hasCandidates(), ((read1->length - read1->Alignments[scoreId2].QStart - read1->Alignments[scoreId2].QEnd)), minResidues, read1->Alignments[scoreId2].Identity >= minIdentity, minIdentity);
		//Log.Message("%d %d %d", read1->length, read1->Alignments[scoreId2].QStart, read1->Alignments[scoreId2].QEnd);
		//Log.Message("Output paired 2: hC %d, R: %d %d, I: %f %f", read2->hasCandidates(), ((read2->length - read2->Alignments[scoreId2].QStart - read2->Alignments[scoreId2].QEnd)), minResidues, read2->Alignments[scoreId2].Identity >= minIdentity, minIdentity);
		//Log.Message("%d %d %d", read2->length, read2->Alignments[scoreId2].QStart, read2->Alignments[scoreId2].QEnd);

		DoWritePair(read1, scoreId1, read2, scoreId2);

		//NGM.AddWrittenRead(read1->ReadId);
		//NGM.AddWrittenRead(read2->ReadId);

	}

	void WriteEpilog() {
		DoWriteEpilog();
	}
};

#endif
