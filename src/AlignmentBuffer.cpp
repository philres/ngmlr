#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>

#include "OutputReadBuffer.h"
#include "Timing.h"

ulong AlignmentBuffer::alignmentCount = 0;
bool AlignmentBuffer::first = true;

void AlignmentBuffer::flush() {
	DoRun();
	nReads = 0;
}

void AlignmentBuffer::debugAlgnFinished(MappedRead * read) {
	Log.Debug(32, "READ_%d\tALGN\tAll alignments computed (%d)", read->ReadId, read->numScores());

	if(read->numScores() > 0) {
		for(int i = 0; i < read->numScores(); ++i) {

			LocationScore score = read->Scores[i];
			Align align = read->Alignments[i];

			SequenceLocation loc = score.Location;
			SequenceProvider.convert(loc);

			int refNameLength = 0;
			//TODO_GENOMESIZE: Re-enable me
			//Log.Debug(128, "READ_%d\tALGN_RESULTS\tCMR_%d\t%f\t%f\t%d\t%s\t%s\t%d\t%s", read->ReadId, i, score.Score.f, align.Identity, align.NM, align.pBuffer1, align.pBuffer2, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength));
		}

	}

#ifdef _DEBUGCMRS
	SequenceLocation rloc = SequenceProvider.convert(cur_read, cur_read->Scores[scoreId].Location.m_Location);
	int refNameLength = 0;
	fprintf(cmrBed, "%s\t%d\t%d\t%s_%d\t%f\t%c\n", SequenceProvider.GetRefName(rloc.getrefId(), refNameLength), rloc.m_Location - (corridor >> 1), rloc.m_Location - (corridor >> 1) + refMaxLen, cur_read->name, scoreId, cur_read->Scores[scoreId].Score.f, (rloc.isReverse()) ? '-' : '+');
#endif
}

void AlignmentBuffer::addRead(MappedRead * read, int scoreID) {
	if (argos) {
		SaveRead(read, read->hasCandidates());
	} else {
		if (!read->hasCandidates() || read->mappingQlty < min_mq) {
			//If read has no CMRs or mapping quality is lower than min mapping quality, output unmapped read
			//read->clearScores(-1);
			SaveRead(read, false);
		} else {
			Log.Debug(512, "READ_%d\tALGN_BUFFER\tCMR_%d %f (location %llu) added to alignment buffer at position %d", read->ReadId, scoreID, read->Scores[scoreID].Score.f, GET_ULOC(read->Scores[scoreID].Location.m_Location), nReads);
			//add alignment computations to buffer. if buffer is full, submit to CPU/GPU
			reads[nReads].scoreId = scoreID;
			reads[nReads++].read = read;
			if (nReads == batchSize) {
				DoRun();
				nReads = 0;
			}
		}
	}
}

void AlignmentBuffer::DoRun() {

	int count = nReads;

	if (count > 0) {
		Log.Debug(32, "INFO\tALGN\tSubmitting %d alignment computations.", count);
		Timer tmr;
		tmr.ST();
		alignmentCount += count;
		for (int i = 0; i < count; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;

			assert(cur_read->hasCandidates());

			//Initialize
			if (cur_read->Scores[scoreID].Location.isReverse()) {
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
			SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
					cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen);

			//initialize arrays for CIGAR and MD string
			static int const qryMaxLen = Config.GetInt("qry_max_len");
			alignBuffer[i].pBuffer1 = new char[std::max(1, qryMaxLen) * 4];
			alignBuffer[i].pBuffer2 = new char[std::max(1, qryMaxLen) * 4];
			*(int*) alignBuffer[i].pBuffer1 = 0x212121;
			*(int*) alignBuffer[i].pBuffer2 = 0x212121;

		}

		//start alignment
		int aligned = aligner->BatchAlign(alignmode | (std::max(outputformat, 1) << 8), count, refBuffer, qryBuffer, alignBuffer,
				(m_EnableBS) ? m_DirBuffer : 0);

		Log.Debug(32, "INFO\tALGN\t%d alignments computed (out of %d)", aligned, count);

		if (aligned != count)
		Log.Error("Error aligning outputs (%i of %i aligned)", aligned, count);

		//process results
		for (int i = 0; i < aligned; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;
			int id = cur_read->ReadId;

			assert(cur_read->hasCandidates());
			cur_read->Scores[scoreID].Location.m_Location += alignBuffer[i].PositionOffset - (corridor >> 1);

			cur_read->Alignments[scoreID] = alignBuffer[i];

			Log.Debug(2048, "READ_%d\tALGN_DETAILS\tCMR_%d\t%f\t%f\t%llu\t%.*s\t%s", cur_read->ReadId, scoreID, cur_read->Scores[scoreID].Score.f, alignBuffer[i].Identity, alignBuffer[i].NM, GET_ULOC(refMaxLen), refBuffer[i], qryBuffer[i]);

			if ((cur_read->Calculated - 1) == scoreID) {

				debugAlgnFinished(cur_read);

				SaveRead(cur_read);
			}

		}
		alignTime = tmr.ET();
	} else {
		Log.Debug(1, "INFO\tALGN\tEmpty buffer submitted.");
	}
}

void AlignmentBuffer::SaveRead(MappedRead * read, bool mapped) {
	//if (!argos) {
		WriteRead(read, mapped);
//	} else {
//		if (mapped) {
//			//Convert mapping position to RefId and position
//			for (int i = 0; i < read->Calculated; ++i) {
//				//TODO: fix for -n > 1
//				//Instead of setting mapped to false set score to 0 and don't print it in the end
//				mapped = SequenceProvider.convert(read->Scores[i].Location);
//			}
//		}
//		OutputReadBuffer::getInstance().addRead(read, mapped);
//		OutputReadBuffer::getInstance().getNextRead(m_Writer);
//	}
}

void AlignmentBuffer::WriteRead(MappedRead* read, bool mapped) {
	static int const topn = Config.GetInt("topn");
	if (mapped) {
		//Convert mapping position to RefId and position
		for (int i = 0; i < read->Calculated; ++i) {
			//TODO: fix for -n > 1
			//Instead of setting mapped to false set score to 0 and don't print it in the end
			mapped = SequenceProvider.convert(read->Scores[i].Location);
		}
	}
	if (read->Paired != 0) {
		if (topn == 1) {
			if (read->Paired->HasFlag(NGMNames::DeletionPending)) {
				if (read->hasCandidates() && read->Paired->hasCandidates()) {
					LocationScore * ls1 = &read->Scores[0];
					LocationScore * ls2 = &read->Paired->Scores[0];
					int distance = ULOC_TO_INT32(
							(ls2->Location.m_Location > ls1->Location.m_Location) ?
									ls2->Location.m_Location - ls1->Location.m_Location + read->length :
									ls1->Location.m_Location - ls2->Location.m_Location + read->Paired->length );

					//int distance = abs(read->TLS()->Location.m_Location - read->Paired->TLS()->Location.m_Location);

					pairInsertCount += 1;
					if (ls1->Location.getrefId() != ls2->Location.getrefId() || distance < _NGM::sPairMinDistance
							|| distance > _NGM::sPairMaxDistance || ls1->Location.isReverse() == ls2->Location.isReverse()) {
						//						Log.Message("%d != %d || %d < _%d || %d > %d || %d == %d", ls1->Location.getrefId() , ls2->Location.getrefId(), distance, _NGM::sPairMinDistance, distance, _NGM::sPairMaxDistance, ls1->Location.isReverse(), ls2->Location.isReverse());
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
