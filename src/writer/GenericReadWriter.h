#ifndef __GENERICREADWRITER_H__
#define __GENERICREADWRITER_H__

#include <map>

#include <stdio.h>
#include <stdarg.h>

#include "ILog.h"
#include "Config.h"
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

		identity = 0;
		if (Config.Exists("identity_tresh")) {
			identity = Config.GetFloat("identity_tresh", 0, 100);
			identity = identity / 100.0f;
		}
		writeUnmapped = !Config.GetInt("no_unal");
	}
	virtual ~GenericReadWriter() {
		if(writeBuffer != 0) {
			delete[] writeBuffer;
			writeBuffer = 0;
		}
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
		va_list arg;

		va_start(arg, format);
		done = vsprintf(writeBuffer + bufferPosition, format, arg);
		bufferPosition += done;
		va_end(arg);
		return done;
	}

private:

public:
	void WriteProlog() {
		DoWriteProlog();
	}

	void WriteRead(MappedRead const * const read, bool mapped = true) {
		if (Config.Exists(ARGOS)) {
			if (mapped) {
				DoWriteRead(read, 0);
				NGM.AddMappedRead(read->ReadId);
			} else {
				DoWriteUnmappedRead(read);
			}
		} else {
			if (mapped) {
				std::map<SequenceLocation, bool> iTable;
				bool mappedOnce = false;
				for (int i = 0; i < read->Calculated; ++i) {

					float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
					float minResidues = Config.GetFloat("min_residues", 0, 1000);

					if (minResidues <= 1.0f) {
						minResidues = read->length * minResidues;
					}

					mapped = mapped && (read->Alignments[i].Identity >= minIdentity);
					mapped = mapped && ((float)(read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd) >= minResidues);

					Log.Debug(4, "READ_%d\tOUTPUT\tChecking alignment CRM_%d (read length: %d - %d - %d)\t%f >= %f\t%f >= %f", read->ReadId, i, read->length, read->Alignments[i].QStart, read->Alignments[i].QEnd, read->Alignments[i].Identity, minIdentity, (float)(read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd), minResidues);

					if (mapped) {
						mappedOnce = true;
						if (iTable.find(read->Scores[i].Location) == iTable.end()) {
							iTable[read->Scores[i].Location] = true;
							Log.Debug(4, "READ_%d\tOUTPUT\tWriting alignment CRM_%d", read->ReadId, i);
							DoWriteRead(read, i);
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
	}

	void WritePair(MappedRead * const read1, int const scoreId1, MappedRead * const read2, int const scoreId2) {
		if(read1->HasFlag(NGMNames::Empty) || read2->HasFlag(NGMNames::Empty)) {
			Log.Debug(LOG_OUTPUT_DETAILS, "Empty read found in pair: %s/%s. Both reads will be discarded and not written to output.", read1->name, read2->name);
		} else {

			//TODO: fix paired end!!! MULTI MAP!

			static float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
			static float const minResiduesConfig = Config.GetFloat("min_residues", 0, 1000);
			float minResidues = minResiduesConfig;
			static int const min_mq = Config.GetInt(MIN_MQ);


			float minResidues1 = 0.0f;
			float minResidues2 = 0.0f;

			if (minResidues <= 1.0f) {
				minResidues1 = read1->length * minResidues;
				minResidues2 = read2->length * minResidues;
			} else {
				minResidues1 = minResidues2 = minResidues;
			}

			bool mapped1 = read1->hasCandidates() && read1->mappingQlty >= min_mq && (read1->Alignments[scoreId1].Identity >= minIdentity)
			&& ((float)(read1->length - read1->Alignments[scoreId1].QStart - read1->Alignments[scoreId1].QEnd) >= minResidues1);
			bool mapped2 = read2->hasCandidates() && read2->mappingQlty >= min_mq && (read2->Alignments[scoreId2].Identity >= minIdentity)
			&& ((float)(read2->length - read2->Alignments[scoreId2].QStart - read2->Alignments[scoreId2].QEnd) >= minResidues2);

			if (!mapped1) {
				read1->clearScores();
			} else {
				NGM.AddMappedRead(read1->ReadId);
			}
			if (!mapped2) {
				read2->clearScores();
			} else {
				NGM.AddMappedRead(read2->ReadId);

			}

			//Log.Message("Output paired 1: hC %d, R: %d %d, I: %f %f", read1->hasCandidates(), ((read1->length - read1->Alignments[scoreId2].QStart - read1->Alignments[scoreId2].QEnd)), minResidues, read1->Alignments[scoreId2].Identity >= minIdentity, minIdentity);
			//Log.Message("%d %d %d", read1->length, read1->Alignments[scoreId2].QStart, read1->Alignments[scoreId2].QEnd);
			//Log.Message("Output paired 2: hC %d, R: %d %d, I: %f %f", read2->hasCandidates(), ((read2->length - read2->Alignments[scoreId2].QStart - read2->Alignments[scoreId2].QEnd)), minResidues, read2->Alignments[scoreId2].Identity >= minIdentity, minIdentity);
			//Log.Message("%d %d %d", read2->length, read2->Alignments[scoreId2].QStart, read2->Alignments[scoreId2].QEnd);

			DoWritePair(read1, scoreId1, read2, scoreId2);
		}
	}

	void WriteEpilog() {
		DoWriteEpilog();
	}
}
;

#endif
