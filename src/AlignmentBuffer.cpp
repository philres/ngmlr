#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>

#include "Timing.h"
#include "Debug.h"

ulong AlignmentBuffer::alignmentCount = 0;
bool AlignmentBuffer::first = true;


void AlignmentBuffer::flush() {
	DoRun();
	nReads = 0;
}

void AlignmentBuffer::addRead(MappedRead * read, int scoreID) {
//	//TODO: remove
//	char const * debugRead = "adb-100bp-20mio-paired.000000558.2";
//	if (strcmp(read->name, debugRead) == 0) {
//		Log.Error("add alignment computation");
//		getchar();
//	}
	if (!read->hasCandidates()) {
		SaveRead(read, false);
	} else {
		reads[nReads].scoreId = scoreID;
		reads[nReads++].read = read;
		if (nReads == batchSize) {
			DoRun();
			nReads = 0;
		}
	}
}

void AlignmentBuffer::DoRun() {

	int count = nReads;
	Timer tmr;
	tmr.ST();

	if (count > 0) {
		Timer tmr;
		tmr.ST();
		alignmentCount += count;
		for (int i = 0; i < count; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;

			if (cur_read->hasCandidates()) {

				// Bei Mapping am Minus-Strang, Position bez�glich +-Strang reporten
				if (cur_read->Scores[scoreID].Location.m_Reverse) {
					qryBuffer[i] = cur_read->RevSeq;

					// RefId auf +-Strang setzen
//					--cur_read->Scores[scoreID].Location.m_RefId;
					if (cur_read->Paired != 0) {
						m_DirBuffer[i] = !(cur_read->ReadId & 1);
					} else {
						m_DirBuffer[i] = 1;
					}

				} else {
					qryBuffer[i] = cur_read->Seq;
					if (cur_read->Paired != 0) {
						m_DirBuffer[i] = cur_read->ReadId & 1; //0 if first pair
					} else {
						m_DirBuffer[i] = 0;
					}
				}

				SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), cur_read->Scores[scoreID].Location.m_RefId,
						cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen);

				static int const qryMaxLen = Config.GetInt("qry_max_len");
				alignBuffer[i].pBuffer1 = new char[std::max(1, qryMaxLen) * 4];
				alignBuffer[i].pBuffer2 = new char[std::max(1, qryMaxLen) * 4];
				*(int*) alignBuffer[i].pBuffer1 = 0x212121;
				*(int*) alignBuffer[i].pBuffer2 = 0x212121;

			} else {
				Log.Warning("Unmapped read submitted to alignment computation!");

				Fatal();
				char * refDummy = const_cast<char *>(refBuffer[i]);
				memset(refDummy, '\0', refMaxLen);
				qryBuffer[i] = dummy; //SequenceProvider.GetQrySequence(cur_read->ReadId);

				alignBuffer[i].pBuffer1 = dBuffer;
				alignBuffer[i].pBuffer2 = dBuffer + dbLen / 2;

			}
			//Log.Message("%s: %s", cur_read->name, refBuffer[i]);
			//Log.Message("%s: %s", cur_read->name, qryBuffer[i]);
		}

		Log.Verbose("Thread %i invoking alignment (count = %i)", 0, count);
		Timer x;
		x.ST();
		int aligned = 0;
//		NGM.AquireOutputLock();
		aligned = aligner->BatchAlign(alignmode | (std::max(outputformat, 1) << 8), count, refBuffer, qryBuffer, alignBuffer,
				(m_EnableBS) ? m_DirBuffer : 0);
//		NGM.ReleaseOutputLock();

		if (aligned == count) {
			Log.Verbose("Output Thread %i finished batch (Size = %i, Elapsed: %.2fs)", 0, count, x.ET());
		} else {
			Log.Error("Error aligning outputs (%i of %i aligned)", aligned, count);
		}

		for (int i = 0; i < aligned; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;
//
//			//TODO: remove
//			char const * debugRead = "adb-100bp-20mio-paired.000000558.2";
//			if (strcmp(cur_read->name, debugRead) == 0) {
//				Log.Error("Alignment computed");
//				getchar();
//			}

			int id = cur_read->ReadId;

			if (cur_read->hasCandidates()) {
				cur_read->Scores[scoreID].Location.m_Location += alignBuffer[i].PositionOffset - (corridor >> 1);
				// TODO: Align liefert keine Scores
//					Log.Message("%s: %d %d %f -> %s, %s", cur_read->name, cur_read->QStart, cur_read->QEnd, cur_read->Identity, cur_read->Buffer1, cur_read->Buffer2);

#ifdef _DEBUGOUT
				Log.Message("Read:   %s", cur_read->name);
				Log.Message("Score: %f", cur_read->Scores[scoreID].Score.f);
				Log.Message("Seq:    %s", cur_read->Seq);
				Log.Message("CIGAR:  %s", cur_read->Buffer1);
				Log.Message("MD:     %s", cur_read->Buffer2);
#endif

				cur_read->Alignments[scoreID] = alignBuffer[i];

				if ((cur_read->Calculated - 1) == scoreID) {
					Log.Verbose("Process aligned read %i,%i (%s)",cur_read->EqualScoringCount, cur_read->ReadId, cur_read->name);
					SaveRead(cur_read);
					NGM.GetReadProvider()->DisposeRead(cur_read);
				}

			} else {
				Log.Error("Unmapped read detected during alignment computation!");
				Fatal();
			}
		} // for
		Log.Verbose("Output Thread %i finished batch in %.2fs", 0, tmr.ET());
		alignTime = tmr.ET();
	} // if count > 0
	else {
		Log.Message("Nothing to do...waiting");
	}
}
