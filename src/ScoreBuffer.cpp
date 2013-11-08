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

void ScoreBuffer::DoRun() {

	if (iScores != 0) {

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
				Log.Warning("Could not decode reference for alignment (read: %s)", cur_read->name);
				Log.Warning("Read sequence: %s", cur_read->Seq);
				memset(const_cast<char *>(m_RefBuffer[i]), 'N', refMaxLen);
			}

			m_ScoreBuffer[i] = -1;

		}

		//Compute scores
		ScoreBuffer::scoreCount += iScores;
		int res = aligner->BatchScore(m_AlignMode, iScores, m_RefBuffer, m_QryBuffer, m_ScoreBuffer, (m_EnableBS) ? m_DirBuffer : 0);
		if (res != iScores)
			Log.Error("SW Kernel couldn't calculate all scores (%i out of %i)", res, iScores);

			//Process results
		brokenPairs = 0;
		for (int i = 0; i < iScores; ++i) {
			MappedRead * cur_read = scores[i].read;
			int scoreId = scores[i].scoreId;


			cur_read->Scores[scoreId].Score.f = m_ScoreBuffer[i];
//			Log.Message("SCORE: %f", cur_read->Scores[scoreId].Score.f);

#ifdef _DEBUGSW
			MappedRead * cur_read = scores[i]->Read;
			SequenceLocation loc = scores[i]->Location;
			SequenceLocation rloc = SequenceProvider.convert(cur_read, loc.m_Location);
			int refNameLength = 0;
			Log.Message("%s - Loc: %u (+), Location: %u (Ref: %s), Score: %f", cur_read->name, loc.m_Location, rloc.m_Location, SequenceProvider.GetRefName(rloc.m_RefId, refNameLength), m_ScoreBuffer[i]);
			//Log.Message("%u %u %u %u", loc.m_Location, corridor, (corridor >> 1), loc.m_Location - (corridor >> 1));
			Log.Message("Strand: %c", (loc.m_RefId & 1) ? '-' : '+');
			Log.Message("Ref:  %.*s", refMaxLen, m_RefBuffer[i]);
			Log.Message("Read: %s", m_QryBuffer[i]);
			getchar();
#endif

			if (!isPaired) {
				if (++cur_read->Calculated == cur_read->numScores()) {
					//all scores computed for single end read
					assert(cur_read->hasCandidates());
					if (topn == 1) {
						top1SE(cur_read);
					} else {
						topNSE(cur_read);
					}
				}
			} else {
				if (++cur_read->Calculated == cur_read->numScores() && cur_read->Paired->Calculated == cur_read->Paired->numScores()) {
					if(strcmp(cur_read->name, "HWUSI-EAS475:1:12:17529:18194#0/1") == 0) {
						//if(strcmp(cur_read->name, "HWUSI-EAS475:1:12:17529:18194#0/1") == 0) {
						Log.Message("FOUND: %s", cur_read->name);
						//Fatal();
					}
					//all scores computed for both mates
					if (topn == 1) {
						if (!fastPairing) {
							if (cur_read->Paired->hasCandidates())
								top1PE(cur_read);
							else
								top1SE(cur_read);
						} else {
							top1SE(cur_read);
							if (cur_read->Paired->hasCandidates())
								top1SE(cur_read->Paired);
						}
					} else {
						topNPE(cur_read);
					}
				}
			}
		}
		scoreTime += tmr.ET();
	}
}

bool sortLocationScore(LocationScore a, LocationScore b) {
	return a.Score.f > b.Score.f;
}

void ScoreBuffer::top1SE(MappedRead* read) {
	static const float minF = -1000000;
	float score_max = minF;
	float score_smax = minF;
	int score_max_loc = 0;
	int score_max_count = 0;
	//count number of top-scoring regions
	for (int j = 0; j < read->numScores(); ++j) {
		if (read->Scores[j].Score.f > score_smax) {
			if (read->Scores[j].Score.f > score_max) {
				score_smax = score_max;
				score_max = read->Scores[j].Score.f;
				score_max_loc = j;
				score_max_count = 1;
			} else if (read->Scores[j].Score.f == score_max) {
				++score_max_count;
				score_smax = score_max;
			} else {
				score_smax = read->Scores[j].Score.f;
			}
		} else if (read->Scores[j].Score.f == score_max) {
			++score_max_count;
		}
	}

	int mq = ceil(MAX_MQ * (score_max - score_smax) / score_max);
	read->mappingQlty = mq;

	if (score_max_count == 1 || !equalOnly) {
		assert(read->hasCandidates());

		//clear all sub-optimal scores and submit best-scoring region to alignment computation
		read->clearScores(score_max_loc);
		read->Calculated = 1;
		read->Alignments = new Align[read->Calculated];
		read->EqualScoringCount = score_max_count;

		out->addRead(read, 0);
	} else {
		//To many equal scoring positions
		read->Calculated = 0;
		read->mappingQlty = 0;
		read->clearScores();

		out->addRead(read, -1);
	}
}

void ScoreBuffer::computeMQ(MappedRead* read) {
	//compute mapping quality
	int mq = MAX_MQ;
	if (read->numScores() > 1) {
		mq = ceil(MAX_MQ * (read->Scores[0].Score.f - read->Scores[1].Score.f) / read->Scores[0].Score.f);
	}
	read->mappingQlty = mq;
}

void ScoreBuffer::topNSE(MappedRead* read) {
	//Sort scores
	std::sort(read->Scores, read->Scores + read->numScores(), sortLocationScore);
	int n = std::min(read->numScores(), topn);
	int j = 1;
	int i = 1;
	//Count number of top-scoring regions
	while (i < n && read->Scores[0].Score.f == read->Scores[i].Score.f) {
		i += 1;
	}
	read->EqualScoringCount = i;
	if (read->EqualScoringCount < topn - 1 || !equalOnly) {
		if (equalOnly) {
			n = i;
		}

		//compute mapping quality
		computeMQ(read);

		read->Calculated = n;
		read->Alignments = new Align[read->Calculated];

		//Submit reads to alignment computation
		out->addRead(read, 0);
		for (; j < n; ++j) {
			if (equalOnly && read->Scores[0].Score.f != read->Scores[j].Score.f) {
				Log.Error("not equal");
				Fatal();
			}
			out->addRead(read, j);
		}
	} else {
		//To many equal scoring positions
		read->mappingQlty = 0;
		read->clearScores();
		out->addRead(read, -1);
	}
}

struct PairScore {
	float score;
	int insertSize;
	int iRead;
	int iMate;
};

void ScoreBuffer::top1PE(MappedRead* read) {

	MappedRead * mate = read->Paired;

	//Sort scores
	std::sort(read->Scores, read->Scores + read->numScores(), sortLocationScore);
	std::sort(mate->Scores, mate->Scores + mate->numScores(), sortLocationScore);

	computeMQ(read);
	computeMQ(mate);

	//Use only scores that are > topScore * cutoff
	float const minScoreRead = read->Scores[0].Score.f * pairScoreCutoff;
	int nScoreRead = 1;
	while (nScoreRead < read->numScores() && minScoreRead <= read->Scores[nScoreRead].Score.f) {
		nScoreRead += 1;
	}

	float const minScoreMate = mate->Scores[0].Score.f * pairScoreCutoff;
	int nScoreMate = 1;
	while (nScoreMate < mate->numScores() && minScoreMate <= mate->Scores[nScoreMate].Score.f) {
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
			if (CheckPairs(&read->Scores[i], read->length, &mate->Scores[j], mate->length, topScore, distance, equalScoreFound)) {
				TopScore1 = i;
				TopScore2 = j;
			}
		}
	}
	Log.Verbose("Pairs with equal score: %d", equalScoreFound);
	if (topScore > 0.0f) {
		if (equalScoreFound <= 0 || !equalOnly) {
			pairDistSum += distance;
			pairDistCount += 1;
			NGM.Stats->insertSize = pairDistSum * 1.0f / (pairDistCount);

			//submit reads to alignment computation
			read->EqualScoringCount = equalScoreFound;
			read->Calculated = 1;
			read->clearScores(TopScore1);
			read->Alignments = new Align[read->Calculated];

			mate->EqualScoringCount = equalScoreFound;
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

bool ScoreBuffer::CheckPairs(LocationScore * ls1, int const readLength1, LocationScore * ls2, int const readLength2, float & pairTopScore, int & insertSize, int & equalScore) {

	//compute insert size
	int currentInsertsize =
			(ls2->Location.m_Location > ls1->Location.m_Location) ?
					ls2->Location.m_Location - ls1->Location.m_Location + readLength2 :
					ls1->Location.m_Location - ls2->Location.m_Location + readLength1;

	if (currentInsertsize > _NGM::sPairMinDistance && currentInsertsize < _NGM::sPairMaxDistance) {
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
		//if(strcmp("HWUSI-EAS475:1:12:17529:18194#0/1", read->name) == 0) {
		//	Log.Error("Read %s: scorebuffer begin %d/%d", read->name, iScores, swBatchSize);
			//getchar();
		//}
		//if(strcmp(read->name, "adb-100bp-20mio-paired.000000071.1") == 0) {
			Log.Verbose("%s\t%u\t%f\tEND", read->name, newScores[i].Location.m_Location, newScores[i].Score.f);
		//}
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

