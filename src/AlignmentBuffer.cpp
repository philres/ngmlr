#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>

#include "Timing.h"

ulong AlignmentBuffer::alignmentCount = 0;
bool AlignmentBuffer::first = true;

void AlignmentBuffer::flush() {
	DoRun();
	nReads = 0;
}

void AlignmentBuffer::addRead(MappedRead * read, int scoreID) {
	if (!read->hasCandidates()) {
		//If read has no CMRs, output unmapped read
		SaveRead(read, false);
	} else {
		//add alignment computations to buffer. if buffer is full, submit to CPU/GPU
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

			assert(cur_read->hasCandidates());

			//Initialize
			if (cur_read->Scores[scoreID].Location.m_Reverse) {
				qryBuffer[i] = cur_read->RevSeq;

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

			//decode reference sequence
			SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), cur_read->Scores[scoreID].Location.m_RefId,
					cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen);

			//initialize arrays for CIGAR and MD string
			static int const qryMaxLen = Config.GetInt("qry_max_len");
			alignBuffer[i].pBuffer1 = new char[std::max(1, qryMaxLen) * 4];
			alignBuffer[i].pBuffer2 = new char[std::max(1, qryMaxLen) * 4];
			*(int*) alignBuffer[i].pBuffer1 = 0x212121;
			*(int*) alignBuffer[i].pBuffer2 = 0x212121;

		}

		Log.Verbose("Thread %i invoking alignment (count = %i)", 0, count);
		//start alignment
		Timer x;
		x.ST();
		int aligned = aligner->BatchAlign(alignmode | (std::max(outputformat, 1) << 8), count, refBuffer, qryBuffer, alignBuffer,
				(m_EnableBS) ? m_DirBuffer : 0);

		if (aligned == count) {
			Log.Verbose("Output Thread %i finished batch (Size = %i, Elapsed: %.2fs)", 0, count, x.ET());
		} else {
			Log.Error("Error aligning outputs (%i of %i aligned)", aligned, count);
		}

		//process results
		for (int i = 0; i < aligned; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;
			int id = cur_read->ReadId;

			assert(cur_read->hasCandidates());
			cur_read->Scores[scoreID].Location.m_Location += alignBuffer[i].PositionOffset - (corridor >> 1);
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
				Log.Verbose("Process aligned read. Equal: %i, ReadId: %i, numScore: %d, calculated: %d, (%s)",cur_read->EqualScoringCount, cur_read->ReadId, cur_read->numScores(), cur_read->Calculated, cur_read->name);
				SaveRead(cur_read);
			}

		}
		Log.Verbose("Output Thread %i finished batch in %.2fs", 0, tmr.ET());
		alignTime = tmr.ET();
	} else {
		Log.Message("Nothing to do...waiting");
	}
}

void AlignmentBuffer::SaveRead(MappedRead* read, bool mapped) {
	static int const topn = Config.GetInt("topn");
	if (mapped) {
		for (int i = 0; i < read->Calculated; ++i) {
			//Get starting positions on the concatenated reference for all chromosomes
			static int refCount = SequenceProvider.GetRefCount();
			if (refStartPos == 0) {
				refStartPos = new int[refCount / ((NGM.DualStrand()) ? 2 : 1)];
				int i = 0;
				int j = 0;
				while (i < refCount/* && loc.m_Location >= SequenceProvider.GetRefStart(i)*/) {
					refStartPos[j++] = SequenceProvider.GetRefStart(i);
					//Log.Message("refstar %d: %d", j-1, SequenceProvider.GetRefStart(i));
					i += (NGM.DualStrand()) ? 2 : 1;
				}
			}

			//Convert position back to Chromosome+Position
					SequenceLocation loc = read->Scores[i].Location;
					//Log.Message("Loc: %u", loc.m_Location);
					int * upper = std::upper_bound(refStartPos, refStartPos + (refCount / ((NGM.DualStrand()) ? 2 : 1)), loc.m_Location);
					//Log.Message("upper %d %d", *upper, *(upper-1));
					std::ptrdiff_t refId = ((upper - 1) - refStartPos) * ((NGM.DualStrand()) ? 2 : 1);
					loc.m_Location -= *(upper - 1);
					loc.m_RefId = refId;
					//Log.Message("Converted score %d: %hd %d %u", i, loc.m_RefId, refId, loc.m_Location);
					read->Scores[i].Location = loc;

					if (loc.m_Reverse) {
						if (read->qlty != 0)
						std::reverse(read->qlty, read->qlty + strlen(read->qlty));
					}
				}
			}
	if (read->Paired != 0) {
		if (topn == 1) {
			if (read->Paired->HasFlag(NGMNames::DeletionPending)) {
				if (read->hasCandidates() && read->Paired->hasCandidates()) {
					LocationScore * ls1 = &read->Scores[0];
					LocationScore * ls2 = &read->Paired->Scores[0];
					int distance =
							(ls2->Location.m_Location > ls1->Location.m_Location) ?
									ls2->Location.m_Location - ls1->Location.m_Location + read->length :
									ls1->Location.m_Location - ls2->Location.m_Location + read->Paired->length;

					//int distance = abs(read->TLS()->Location.m_Location - read->Paired->TLS()->Location.m_Location);

					pairInsertCount += 1;
					if (ls1->Location.m_RefId != ls2->Location.m_RefId || distance < _NGM::sPairMinDistance
							|| distance > _NGM::sPairMaxDistance || ls1->Location.m_Reverse == ls2->Location.m_Reverse) {
						read->SetFlag(NGMNames::PairedFail);
						read->Paired->SetFlag(NGMNames::PairedFail);
						brokenPairs += 1;
					} else {
						pairInsertSum += distance;
					}
				}
				m_Writer->WritePair(read, 0, read->Paired, 0);
			}
		} else {
			Log.Error("TopN > 1 is currently not supported for paired end reads.");
			Fatal();
		}
	} else {
		m_Writer->WriteRead(read, mapped);
	}
	if (pairInsertCount % 1000 == 0) {
		NGM.Stats->validPairs = (pairInsertCount - brokenPairs) * 100.0f / pairInsertCount;
//		NGM.Stats->insertSize = tSum * 1.0f / (tCount - brokenPairs);
	}
	NGM.GetReadProvider()->DisposeRead(read);
}
