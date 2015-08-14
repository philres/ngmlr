/*
 * SWwoBuffer.cpp
 *
 *  Created on: Feb 19, 2013
 *      Author: philipp_
 */

#include "ScoreBuffer.h"

#include <algorithm>
#include <cmath>

#include "NGMTask.h"
#include "Debug.h"
#include "Timing.h"
#include "CS.h"
#include "AlignmentBuffer.h"

ulong ScoreBuffer::scoreCount = 0;

float const MAX_MQ = 60.0f;

struct PairScore {
	float score;
	int insertSize;
	int iRead;
	int iMate;
};

bool sortLocationScore(LocationScore a, LocationScore b) {
	return a.Score.f > b.Score.f;
}

int ScoreBuffer::computeMQ(float bestScore, float secondBestScore) {
	int mq = ceil(MAX_MQ * (bestScore - secondBestScore) / bestScore);
	return mq;
}

void ScoreBuffer::computeMQ(MappedRead* read) {
	//compute mapping quality
	int mq = MAX_MQ;
	if (read->numScores() > 1) {
		mq = computeMQ(read->Scores[0].Score.f, read->Scores[1].Score.f);
	}
	read->mappingQlty = mq;
}

void ScoreBuffer::debugScoresFinished(MappedRead * read) {
	Log.Debug(16, "READ_%d\tSCORES\tAll scores computed (%d)", read->ReadId, read->numScores());

	if(read->numScores() > 0) {
		LocationScore * tmpScores = new LocationScore[read->numScores()];
		memcpy(tmpScores, read->Scores, sizeof(LocationScore) * read->numScores());

		//TODO: sort by location
		std::sort(tmpScores, tmpScores + read->numScores(), sortLocationScore);

		for(int i = 0; i < read->numScores(); ++i) {
			LocationScore score = tmpScores[i];

			SequenceLocation loc = score.Location;
			SequenceProvider.convert(loc);

			int refNameLength = 0;
			Log.Debug(64, "READ_%d\tSCORES_RESULTS\tCMR_%d\t%f\t%llu\t%s", read->ReadId, i, score.Score.f, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength));
		}

	}

#ifdef _DEBUGCMRS
	SequenceLocation rloc = SequenceProvider.convert(cur_read, cur_read->Scores[scoreId].Location.m_Location);
	int refNameLength = 0;
	fprintf(cmrBed, "%s\t%d\t%d\t%s_%d\t%f\t%c\n", SequenceProvider.GetRefName(rloc.getrefId(), refNameLength), rloc.m_Location - (corridor >> 1), rloc.m_Location - (corridor >> 1) + refMaxLen, cur_read->name, scoreId, cur_read->Scores[scoreId].Score.f, (rloc.isReverse()) ? '-' : '+');
#endif
}

ReadGroup* ScoreBuffer::updateGroupInfo(MappedRead* cur_read) {
	ReadGroup* group = cur_read->group;
	//TODO: make atomic, parts from a group can end up in different threads!
	group->readsFinished += 1;
	//				Log.Message("Scorecount for %s (%d): %d - %d", cur_read->name, cur_read->ReadId, cur_read->numScores(), cur_read->Calculated);
	if (cur_read->Scores[0].Location.isReverse()) {
		group->reverseMapped += 1;
	} else {
		group->fwdMapped += 1;
	}
	group->bestScoreSum += (int) (cur_read->Scores[0].Score.f);
	return group;
}

void ScoreBuffer::DoRun() {

	if (iScores != 0) {
		Log.Debug(16, "INFO\tSCORES\tSubmitting %d score computations.", iScores);
		Timer tmr;
		tmr.ST();
		//Prepare for score computation
		for (int i = 0; i < iScores; ++i) {
			MappedRead * cur_read = scores[i].read;
			int scoreId = scores[i].scoreId;

			//Initialize buffers for score computation
			SequenceLocation loc = cur_read->Scores[scoreId].Location;
			if (NGM.DualStrand() && (loc.isReverse())) {
				// RefId auf +-Strang setzen
				//--loc.m_RefId;
				cur_read->computeReverseSeq();
				m_QryBuffer[i] = cur_read->RevSeq;
				if (isPaired) {
					m_DirBuffer[i] = !(cur_read->ReadId & 1);
				} else {
					m_DirBuffer[i] = 1;
				}
			} else {
				m_QryBuffer[i] = cur_read->Seq;
				if (isPaired) {
					m_DirBuffer[i] = cur_read->ReadId & 1; //0 if first pair
				} else {
					m_DirBuffer[i] = 0;
				}
			}

			//decode reference sequence
			if (!SequenceProvider.DecodeRefSequence(const_cast<char *>(m_RefBuffer[i]), 0,
							loc.m_Location - (corridor >> 1), refMaxLen)) {
//							loc.m_Location - (corridor >> 1), cur_read->length + corridor)) {
				Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", loc.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
				//Log.Warning("Read sequence: %s", cur_read->Seq);
				memset(const_cast<char *>(m_RefBuffer[i]), 'N', refMaxLen);
			}

			//Log.Message("Ref: %s\nSeq: %s\n", m_RefBuffer[i], m_QryBuffer[i]);
			m_ScoreBuffer[i] = -1;

		}

		//Compute scores
		//TODO: move to NGMStats
		ScoreBuffer::scoreCount += iScores;
		int res = aligner->BatchScore(m_AlignMode, iScores, m_RefBuffer, m_QryBuffer, m_ScoreBuffer, (m_EnableBS) ? m_DirBuffer : 0);

		Log.Debug(16, "INFO\tSCORES\t%d scores computed (out of %d)", res, iScores);

		if (res != iScores)
		Log.Error("Kernel couldn't calculate all scores (%i out of %i)", res, iScores);

		//Process results
		brokenPairs = 0;
		for (int i = 0; i < iScores; ++i) {

			MappedRead * cur_read = scores[i].read;
			int scoreId = scores[i].scoreId;
			cur_read->Scores[scoreId].Score.f = m_ScoreBuffer[i];

			//TODO_GENOMESIZE: Re-enable me
			//Log.Debug(1024, "READ_%d\tSCORES_DETAILS\tCMR_%d\t%f\t%.*s\t%s", cur_read->ReadId, scoreId, m_ScoreBuffer[i], refMaxLen, m_RefBuffer[i], m_QryBuffer[i]);

			if (++cur_read->Calculated == cur_read->numScores()) {
				//all scores computed for single end read
				assert(cur_read->hasCandidates());

#ifdef DEBUGLOG
				debugScoresFinished(cur_read);
#endif

				topNSE(cur_read);


				ReadGroup* group = updateGroupInfo(cur_read);

				//If all reads from group are finished
				if(group->readsFinished == group->readNumber) {
					out->processLongReadLIS(group);
				}

			}
		}
		scoreTime += tmr.ET();
	} else {
		Log.Debug(16, "INFO\tSCORES\tEmpty buffer submitted.");
	}
}

void ScoreBuffer::top1SE(MappedRead* read) {

	float bestScore = 0.0f;
	float secondBestScore = 0.0f;

	int bestScoreIndex = 0;
	int numBestScore = 0;

//Find top-scoring CMR
	for (int j = 0; j < read->numScores(); ++j) {
		if (read->Scores[j].Score.f > secondBestScore) {
			if (read->Scores[j].Score.f > bestScore) {
				secondBestScore = bestScore;
				bestScore = read->Scores[j].Score.f;
				bestScoreIndex = j;
				numBestScore = 1;
			} else if (read->Scores[j].Score.f == bestScore) {
				++numBestScore;
				secondBestScore = bestScore;
			} else {
				secondBestScore = read->Scores[j].Score.f;
			}
		} else if (read->Scores[j].Score.f == bestScore) {
			++numBestScore;
		}
	}

	read->mappingQlty = computeMQ(bestScore, secondBestScore);

	Log.Debug(16, "READ_%d\tSCORES\tBest score %f (CMR_%d) number %d, second best %f, MQ: %d", read->ReadId, bestScore, bestScoreIndex, numBestScore, secondBestScore, read->mappingQlty);

	if (numBestScore == 1 || !topScoresOnly) {
		assert(read->hasCandidates());

		//clear all sub-optimal scores and submit best-scoring region to alignment computation
		read->clearScores(bestScoreIndex);
		read->Calculated = 1;
		read->Alignments = new Align[read->Calculated];
		read->numTopScores = numBestScore;

		out->addRead(read, 0);
	} else {
		//To many equal scoring positions
		read->Calculated = 0;
		read->mappingQlty = 0;
		read->clearScores();

		out->addRead(read, -1);
	}
}

void ScoreBuffer::topNSE(MappedRead* read) {

//Sort scores
	std::sort(read->Scores, read->Scores + read->numScores(),
			sortLocationScore);

	int numScores = read->numScores();

	int numTopScores = 1;
//Count number of top-scoring regions
	while (numTopScores < numScores
			&& read->Scores[0].Score.f == read->Scores[numTopScores].Score.f) {
		numTopScores += 1;
	}

	read->numTopScores = numTopScores;

	if (read->numTopScores <= maxTopScores || !topScoresOnly) {

		//Set number of scores to number of best scores (only these will be reported)
		if (topScoresOnly) {
			numScores = numTopScores;
		}

		//compute mapping quality
		computeMQ(read);

		//numScores alignments will be computed in the next step
		read->Calculated = numScores;
		read->Alignments = new Align[read->Calculated];

//		//Submit reads to alignment computation
//		out->addRead(read, 0);
//		for (int j = 1; j < numScores; ++j) {
//			if (topScoresOnly
//					&& read->Scores[0].Score.f != read->Scores[j].Score.f) {
//				Log.Error("Internal error while processing alignment scores for read %s", read->name);
//				Fatal();
//			}
//			out->addRead(read, j);
//		}
	} else {
//		//To many equal scoring positions, report as unmapped
//		read->mappingQlty = 0;
//		read->clearScores();
//		out->addRead(read, -1);
	}
}

//void ScoreBuffer::topNSE(MappedRead* read) {
//	//Sort scores
//	std::sort(read->Scores, read->Scores + read->numScores(), sortLocationScore);
//	int n = std::min(read->numScores(), topn);
//	int j = 1;
//	int i = 1;
//	//Count number of top-scoring regions
//	while (i < n && read->Scores[0].Score.f == read->Scores[i].Score.f) {
//		i += 1;
//	}
//	read->EqualScoringCount = i;
//	if (read->EqualScoringCount < topn - 1 || !equalOnly) {
//		if (equalOnly) {
//			n = i;
//		}
//
//		//compute mapping quality
//		computeMQ(read);
//
//		read->Calculated = n;
//		read->Alignments = new Align[read->Calculated];
//
//		//Submit reads to alignment computation
//		out->addRead(read, 0);
//		for (; j < n; ++j) {
//			if (equalOnly && read->Scores[0].Score.f != read->Scores[j].Score.f) {
//				Log.Error("not equal");
//				Fatal();
//			}
//			out->addRead(read, j);
//		}
//	} else {
//		//To many equal scoring positions
//		read->mappingQlty = 0;
//		read->clearScores();
//		out->addRead(read, -1);
//	}
//}

void ScoreBuffer::top1PE(MappedRead* read) {

	MappedRead * mate = read->Paired;

	//Sort scores
	std::sort(read->Scores, read->Scores + read->numScores(),
			sortLocationScore);
	std::sort(mate->Scores, mate->Scores + mate->numScores(),
			sortLocationScore);

	computeMQ(read);
	computeMQ(mate);

	//Use only scores that are > topScore * cutoff
	float const minScoreRead = read->Scores[0].Score.f * pairScoreCutoff;
	int nScoreRead = 1;
	while (nScoreRead < read->numScores()
			&& minScoreRead <= read->Scores[nScoreRead].Score.f) {
		nScoreRead += 1;
	}

	float const minScoreMate = mate->Scores[0].Score.f * pairScoreCutoff;
	int nScoreMate = 1;
	while (nScoreMate < mate->numScores()
			&& minScoreMate <= mate->Scores[nScoreMate].Score.f) {
		nScoreMate += 1;
	}

	//Look for top-scoring pair. If pairs with equal score are found ties are broken by comparing their
	//insert size to the average insert size
	float topScore = 0.0f;

	int distance = 0;
	int equalScoreFound = 0;

	int TopScore1 = -1;
	int TopScore2 = -1;

	Log.Verbose("Read: %d (%d), Mate (%d): %d", nScoreRead, read->numScores(), nScoreMate, mate->numScores());

	for (int i = 0; i < nScoreRead; ++i) {
		for (int j = 0; j < nScoreMate; ++j) {
			if (CheckPairs(&read->Scores[i], read->length, &mate->Scores[j],
					mate->length, topScore, distance, equalScoreFound)) {
				TopScore1 = i;
				TopScore2 = j;
			}
		}
	}
	Log.Verbose("Pairs with equal score: %d", equalScoreFound);
	if (topScore > 0.0f) {
		if (equalScoreFound <= 0 || !topScoresOnly) {
			pairDistSum += distance;
			pairDistCount += 1;
			NGM.Stats->insertSize = pairDistSum * 1.0f / (pairDistCount);

			//submit reads to alignment computation
			read->numTopScores = equalScoreFound;
			read->Calculated = 1;
			read->clearScores(TopScore1);
			read->Alignments = new Align[read->Calculated];

			mate->numTopScores = equalScoreFound;
			mate->Calculated = 1;
			mate->clearScores(TopScore2);
			mate->Alignments = new Align[mate->Calculated];

			out->addRead(read, 0);
			out->addRead(mate, 0);
		} else {
			//To many equal scoring positions
			read->mappingQlty = 0;
			read->clearScores();
			mate->mappingQlty = 0;
			mate->clearScores();

			out->addRead(read, -1);
			out->addRead(mate, -1);
		}

	} else {
		//No proper pair found, map single-end.
		read->SetFlag(NGMNames::PairedFail);
		mate->SetFlag(NGMNames::PairedFail);

		top1SE(read);
		top1SE(mate);
	}
}

void ScoreBuffer::topNPE(MappedRead* read) {
	Log.Error("Paired end mode with topn > 1 not yet supported.");
	Fatal();
}

bool ScoreBuffer::CheckPairs(LocationScore * ls1, int const readLength1,
		LocationScore * ls2, int const readLength2, float & pairTopScore,
		int & insertSize, int & equalScore) {

	//compute insert size
	int currentInsertsize =
			(ls2->Location.m_Location > ls1->Location.m_Location) ?
					ls2->Location.m_Location - ls1->Location.m_Location
							+ readLength2 :
					ls1->Location.m_Location - ls2->Location.m_Location
							+ readLength1;

	if (currentInsertsize > _NGM::sPairMinDistance
			&& currentInsertsize < _NGM::sPairMaxDistance) {
		//Log.Green("[%d, %d] Score: %f, Distance: %d (%u - %u), AVG: %d", ls1->Read->ReadId, ls2->Read->ReadId, ls1->Score.f + ls2->Score.f, distance, ls2->Location.m_Location, ls1->Location.m_Location, tSum / tCount);
		float pairScore = ls1->Score.f + ls2->Score.f;
		if (pairScore > pairTopScore * 1.00f) {
			//New top-scoring pair found
			pairTopScore = pairScore;
			insertSize = currentInsertsize;
			return true;
		} else if (pairScore == pairTopScore /** pairScoreCutoff*/) {
			//pair with equal score found
			Log.Verbose("Found a pair (%d, %d) with same score (%f). Choosing pair according to insert size.", ls1->Read->ReadId, ls2->Read->ReadId, pairScore);
			int avg = pairDistSum / pairDistCount;
			if (abs(insertSize - avg) > abs(currentInsertsize - avg)) {
				//choosing pair that is closer to average insert size
				pairTopScore = pairScore;
				insertSize = currentInsertsize;
				return true;
			} else if (abs(insertSize) == abs(currentInsertsize)) {
				//pairs have same insert size, choosing randomly
				equalScore += 1;
				Log.Verbose("Same insert size (%d). Choosing randomly.", currentInsertsize);
			}
		}
	}

	return false;
}

void ScoreBuffer::addRead(MappedRead * read, int count) {

	LocationScore * newScores = read->Scores;
	if (count == 0) {
		Log.Error("Internal error (count == 0). Please report this on https://github.com/Cibiv/NextGenMap/issues");
		Fatal();
	}

	//Adding scores to buffer. If buffer full, submit to CPU/GPU for score computation
	for (int i = 0; i < count; ++i) {

		Log.Debug(256, "READ_%d\tSCORES_BUFFER\tCMR_%d %f (location %llu) added to score buffer at position %d", read->ReadId, i, newScores[i].Score.f, newScores[i].Location.m_Location, iScores);

		scores[iScores].read = read;
		scores[iScores++].scoreId = i;
		if(iScores == swBatchSize) {
			//Log.Error("Read %s: scorebuffer begin %d/%d", read->name, iScores, swBatchSize);
			DoRun();
			iScores = 0;
		}
	}
}

void ScoreBuffer::flush() {
	//Force submitting remaining computation from buffer to CPU/GPU
	DoRun();
	iScores = 0;
	out->flush();
}

