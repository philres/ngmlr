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
			Log.Debug(512, "READ_%d\tALGN_BUFFER\tCMR_%d %f (location %llu) added to alignment buffer at position %d", read->ReadId, scoreID, read->Scores[scoreID].Score.f, read->Scores[scoreID].Location.m_Location, nReads);
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
			if (!SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
							cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen)) {
//							cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor)) {
				Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
				//Log.Warning("Read sequence: %s", cur_read->Seq);
				memset(const_cast<char *>(refBuffer[i]), 'N', refMaxLen);
			}
//			//decode reference sequence
//			SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
//					cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen);

			//initialize arrays for CIGAR and MD string
			static int const qryMaxLen = Config.GetInt("qry_max_len");
			alignBuffer[i].pBuffer1 = new char[std::max(1, qryMaxLen) * 4];
			alignBuffer[i].pBuffer2 = new char[std::max(1, qryMaxLen) * 4];
			*(int*) alignBuffer[i].pBuffer1 = 0x212121;
			*(int*) alignBuffer[i].pBuffer2 = 0x212121;

			//Log.Message("Ref:  %s\nRead: %s", refBuffer[i], qryBuffer[i]);

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

			Log.Debug(2048, "READ_%d\tALGN_DETAILS\tCMR_%d\t%f\t%f\t%llu\t%.*s\t%s", cur_read->ReadId, scoreID, cur_read->Scores[scoreID].Score.f, alignBuffer[i].Identity, alignBuffer[i].NM, refMaxLen, refBuffer[i], qryBuffer[i]);

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

Align AlignmentBuffer::computeAlignment(uloc const position, int const corridor,
		char * const readSeq, size_t const readLength) {

	Align align;

	size_t const refSeqLen = readLength + corridor + 1;
	char * refSeq = new char[refSeqLen];

	//decode reference sequence
	if (!SequenceProvider.DecodeRefSequenceExact(refSeq, position, refSeqLen, corridor)) {
		//Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
		Log.Warning("Could not decode reference for alignment");
		memset(refSeq, 'N', refSeqLen);
	}
	//initialize arrays for CIGAR and MD string
	align.pBuffer1 = new char[readLength * 4];
	align.pBuffer2 = new char[readLength * 4];
	*(int*) align.pBuffer1 = 0x212121;
	*(int*) align.pBuffer2 = 0x212121;

	//Local alignment
	int mode = 0;
	Log.Message("Aligning %d bp to %d bp (corridor %d)", readLength, refSeqLen, corridor);
	printf("%llu\t%d\t%s\t%s", position, corridor, refSeq, readSeq);
	aligner->SingleAlign(mode, corridor, (char const * const ) refSeq,
			(char const * const ) readSeq, align, 0);
//	Log.Message(">Ref\n%s\n>Read\n%s", refSeq, readSeq);

	delete[] refSeq;
	refSeq = 0;

	return align;
}

Align AlignmentBuffer::computeAlignment(MappedRead* read, int const scoreId,
		int const corridor) {

	Align align;
	LocationScore & score = read->Scores[scoreId];

	Log.Message("Computing alignment (%d) for position: %llu", scoreId, score.Location.m_Location);
	Log.Message("Corridor: %d", corridor);
	char * refBuffer = new char[read->length + corridor + 1];

	//decode reference sequence
	if (!SequenceProvider.DecodeRefSequence(refBuffer, 0,
			score.Location.m_Location - (corridor >> 1), read->length + corridor)) {
		//Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
		Log.Warning("Could not decode reference for alignment (read: %s)", read->name);
		memset(refBuffer, 'N', read->length * 1.2f);
	}
	//initialize arrays for CIGAR and MD string
	align.pBuffer1 = new char[read->length * 4];
	align.pBuffer2 = new char[read->length * 4];
	*(int*) align.pBuffer1 = 0x212121;
	*(int*) align.pBuffer2 = 0x212121;

	int const mode = 0;
	if (score.Location.isReverse()) {
		aligner->SingleAlign(mode, corridor, (char const * const ) refBuffer,
				(char const * const ) read->RevSeq, align, 0);
//		printf(">Ref_%s\n%s\n>%s_rev\n%s", read->name, refBuffer, read->name,
//				read->RevSeq);
	} else {
		aligner->SingleAlign(mode, corridor, (char const * const ) refBuffer,
				(char const * const ) read->Seq, align, 0);
//		printf(">Ref_%s\n%s\n>%s\n%s", read->name, refBuffer, read->name,
//				read->Seq);
	}

	score.Location.m_Location += align.PositionOffset - (corridor >> 1);

	delete[] refBuffer;

	return align;
}

void AlignmentBuffer::processLongRead(ReadGroup * group) {
	//					Log.Message("Read group with id %d finished", group->readId);
	//					Log.Message("Name: %s", cur_read->name);
	//					Log.Message("Reads in group: %d", group->readNumber);
	//					Log.Message("Reads finished: %d", group->readsFinished);
	//					Log.Message("Fwd: %d, Rev: %d", group->fwdMapped, group->reverseMapped);
	//					Log.Message("Avg best score %f", group->bestScoreSum * 1.0f / group->readsFinished);

	float avgGroupScore = group->bestScoreSum * 1.0f / group->readsFinished;
	float minGroupScore = avgGroupScore * 0.8f;

//	for (int j = 0; j < group->readNumber; ++j) {
//		MappedRead * part = group->reads[j];
//		//Log.Message("ID: %d (has %d scores)", part->ReadId, part->numScores());
//		float minScore = part->Scores[0].Score.f * 0.8;
//		for (int k = 0; k < part->numScores(); ++k) {
//			if (part->Scores[k].Score.f > minScore) {
//				//Log.Message("\t%f at %llu", part->Scores[k].Score.f, part->Scores[k].Location.m_Location);
//			}
//		}
//	}

	//Find first read part that maps with a min score
	int first = 0;
	while ((group->reads[first]->numScores() == 0
			|| group->reads[first]->Scores[0].Score.f < minGroupScore)
			&& first < group->readNumber) {
		first += 1;
	}

	//Find last read part that maps with a min score
	int last = group->readNumber - 1;
	while ((group->reads[last]->numScores() == 0
			|| group->reads[last]->Scores[0].Score.f < minGroupScore)
			&& last >= 0) {
		last -= 1;
	}

	if (first == group->readNumber || last < 0) {
		Log.Message("Could not map read.");
	} else {
		//Distance on read between start of first and last mapped read part
		//+1 to take the full last read part into account
		int distOnRead = (last - first + 1) * 512;

		//If not the whole read is aligned add half of the read part size to alignment
		size_t startPosOnRead = std::max(0, first * 512 - 256);
		size_t endPosOnRead = std::min(group->fullRead->length, (last + 1) * 512 + 256);

		Log.Message("Name: %s", group->fullRead->name);
		Log.Message("On read (length %d): %d to %d (dist %d)", group->fullRead->length, startPosOnRead, endPosOnRead, distOnRead);

		bool isReverse = false;
		//Compute distance between mapped location of first and last mapped read part on reference
		uloc endPos = group->reads[last]->Scores[0].Location.m_Location;
		uloc startPos = group->reads[first]->Scores[0].Location.m_Location;
		int distOnRef = 0;
		if(startPos > endPos) {
			distOnRef = startPos - endPos;
			isReverse = true;
		} else {
			distOnRef = endPos - startPos;
		}
		//+512 to take full length of last read part into account
		distOnRef += 512;
		Log.Message("Start pos on ref: %llu to %llu (dist %llu)", startPos, endPos, distOnRef);

		//Log.Message("Start: %llu, End: %llu", first, last);
		//Log.Message("Read length: %d", group->fullRead->length);

		float coveredOnRead = (distOnRead) * 100.0f / group->fullRead->length;
		//Log.Message("Covered on read: %f", coveredOnRead);
		//Log.Message("On read: %d, On ref: %llu", distOnRead, distOnRef);

		//Difference between distance in read cooridnates and distance in ref coordinates
		//If read doesn't span larger structural variations, difference should only be
		//caused by PacBio sequence error model
		int difference = distOnRead - distOnRef;
		float diffPerc = (distOnRead - distOnRef) * 1.0f / group->fullRead->length;
		//Log.Message("Difference: %d (%f)", difference, diffPerc);

//		printf("%s\t%d\t%d\%d\t%d\t%d\t%d\t%d\t%f\t%d\t%f\n", group->fullRead->name,
//				group->fullRead->ReadId,
//				group->fullRead->length,
//				group->readNumber,
//				first, last,
//				distOnRead, distOnRef,
//				coveredOnRead,
//				difference, diffPerc);

		//If difference < 0.1. assume that read doesn't span larger SVs and map
		//using alignment
		if(abs(diffPerc) < 0.1) {
			//Normal read. No event.

			MappedRead * read = group->fullRead;

			LocationScore * tmp = new LocationScore();

			if(startPos > endPos) {
				tmp->Location.m_Location = endPos;
				tmp->Location.setReverse(true);
			} else {
				tmp->Location.m_Location = startPos;
				tmp->Location.setReverse(false);
			}

			read->AllocScores(tmp, 1);
			read->Alignments = new Align[1];

			Timer tmr;
			tmr.ST();
			int corridor = std::max(abs(difference) * 2, (int)(read->length * 0.05));

			if(isReverse) {
				Log.Message("Read mapped reverse");

				read->computeReverseSeq();
				char * const readSeq = read->RevSeq + (read->length - endPosOnRead);
				size_t const readSeqLen = endPosOnRead - startPosOnRead;
				Log.Message("ReadSeqLen: %d", readSeqLen);

				printf("%s\t", read->name);
				read->Alignments[0] = computeAlignment(read->Scores[0].Location.m_Location, corridor, readSeq, readSeqLen);
			} else {
				char * const readSeq = read->Seq + startPosOnRead;
				size_t const readSeqLen = endPosOnRead - startPosOnRead;

				printf("%s\t", read->name);
				read->Alignments[0] = computeAlignment(read->Scores[0].Location.m_Location, corridor, readSeq, readSeqLen);
			}
			Log.Message("CIGAR: %s", read->Alignments[0].pBuffer1);

			printf("\t%d\n", read->Alignments[0].PositionOffset);
			read->Scores[0].Location.m_Location += read->Alignments[0].PositionOffset - (corridor >> 1);

			Log.Message("Alignment took %fs", tmr.ET());

			read->Calculated = 1;

			WriteRead(read, true);

		} else {
			WriteRead(group->fullRead, false);
		}
	}
//	getchar();
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
					int distance =
							(ls2->Location.m_Location > ls1->Location.m_Location) ?
									ls2->Location.m_Location
											- ls1->Location.m_Location
											+ read->length :
									ls1->Location.m_Location
											- ls2->Location.m_Location
											+ read->Paired->length;

					//int distance = abs(read->TLS()->Location.m_Location - read->Paired->TLS()->Location.m_Location);

					pairInsertCount += 1;
					if (ls1->Location.getrefId() != ls2->Location.getrefId()
							|| distance < _NGM::sPairMinDistance
							|| distance > _NGM::sPairMaxDistance
							|| ls1->Location.isReverse()
									== ls2->Location.isReverse()) {
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
