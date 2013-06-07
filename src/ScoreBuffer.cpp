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

//static int refMaxLen = 0;
//static int qryMaxLen = 0;

float const MAX_MQ = 60.0f;

long tCount = 1;
long tSum = 0;
long brokenPairs = 0;

void ScoreBuffer::DoRun() {

	if (iScores != 0) {

		Timer tmr;
		tmr.ST();

		for (int i = 0; i < iScores; ++i) {
			MappedRead * cur_read = scores[i]->Read;

//			//TODO: remove
//			char const * debugRead = "adb-100bp-20mio-paired.000000558.2";
//			if (strcmp(cur_read->name, debugRead) == 0) {
//				Log.Error("Queueing score computation.");
//				getchar();
//			}

			SequenceLocation loc = scores[i]->Location;
			if (NGM.DualStrand() && (loc.m_RefId & 1)) {
				// RefId auf +-Strang setzen
				--loc.m_RefId;
				cur_read->computeReverseSeq();
				m_QryBuffer[i] = cur_read->RevSeq;
				if (NGM.Paired()) {
					m_DirBuffer[i] = !(cur_read->ReadId & 1);
				} else {
					m_DirBuffer[i] = 1;
				}
			} else {
				m_QryBuffer[i] = cur_read->Seq;
				if (NGM.Paired()) {
					m_DirBuffer[i] = cur_read->ReadId & 1; //0 if first pair
				} else {
					m_DirBuffer[i] = 0;
				}
			}

			if (!SequenceProvider.DecodeRefSequence(const_cast<char *>(m_RefBuffer[i]), loc.m_RefId,
					loc.m_Location - (corridor >> 1), refMaxLen)) {
				Log.Warning("Could not decode reference for alignment (read: %s)", cur_read->name);
				Log.Warning("Read sequence: %s", cur_read->Seq);
				memset(const_cast<char *>(m_RefBuffer[i]), 'N', refMaxLen);
			}

			m_ScoreBuffer[i] = -1;

		}
//				Log.Message("SW Thread %i launching batch (size=%i)", m_TID, count);
//		Timer tmr;
//		tmr.ST();
		ScoreBuffer::scoreCount += iScores;
		int res = 0;
//		NGM.AquireOutputLock();
		res = aligner->BatchScore(m_AlignMode, iScores, m_RefBuffer, m_QryBuffer, m_ScoreBuffer, (m_EnableBS) ? m_DirBuffer : 0);
//		NGM.ReleaseOutputLock();
		if (res != iScores)
			Log.Error("SW Kernel couldnt calculate all scores (%i out of %i)", res, iScores);

//				Log.Warning("SW Thread %i finished batch (Size = %i, Elapsed: %.2fs)", m_TID, count, tmr.ET());
//		Timer y;
//		y.ST();
		brokenPairs = 0;
		for (int i = 0; i < iScores; ++i) {

//			//TODO: remove
//			char const * debugRead = "adb-100bp-20mio-paired.000000558.2";
//			if (strcmp(scores[i]->Read->name, debugRead) == 0) {
//				Log.Error("Score: %f", m_ScoreBuffer[i]);
//				getchar();
//			}

			scores[i]->Score.f = m_ScoreBuffer[i];

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

			//TODO: remove
//			if (strcmp(scores[i]->Read->name, debugRead) == 0) {
//				Log.Error("%d %d %d", scores[i]->Read->Calculated, ++scores[i]->Read->Calculated, scores[i]->Read->numScores());
//				scores[i]->Read->Calculated -= 1;
//			}

			if (++scores[i]->Read->Calculated == scores[i]->Read->numScores()) {
//			if (AtomicInc(&(scores[i]->Read->Calculated)) == scores[i]->Read->numScores()) {
				SendToPostprocessing(scores[i]->Read);
			}
		}
		scoreTime += tmr.ET();
	}
}

bool sortLocationScore(LocationScore a, LocationScore b) {
	return a.Score.f > b.Score.f;
}

void ScoreBuffer::SendToPostprocessing(MappedRead * read) {

	// Program runs in Paired mode and current read got a pair
	Log.Verbose("[SINGLE] SW::SendToPostprocessing: %i (%s)", read->ReadId, read->name);
	if (read->hasCandidates()) {

		static int const topn = Config.GetInt("topn");

		if(topn == 1) {
			float score_max = 0;
			float score_smax = 0;
			int score_max_loc = 0;
			int score_max_count = 0;

			// zähle anzahl topscores
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

			//static const int skip = (Config.Exists("kmer_skip") ? Config.GetInt("kmer_skip", 0, -1) : 0) + 1;
			//float max = (read->length - CS::prefixBasecount + 1) / skip;

			//float t = 100.0f * (read->s / max);
			//float t = 100.0f * (score_max - score_smax) / score_max;
			//int mq = ceil(log(t * 10.0f + 1) / log(11.0f) * 60.0f);
			//int mq = ceil(100.0f * (read->s / max));

			int mq = ceil(MAX_MQ * (score_max - score_smax) / score_max);
//		Log.Message("%s: %f %f -> %d, (%d) => %d", read->name, score_max, score_smax, mq, read->mappingQlty, mq2);
			//read->mappingQlty = std::min(mq, read->mappingQlty);
			read->mappingQlty = mq;

			//TODO: fix SAM tag X0
			//read->EqualScoringCount = score_max_count;
			read->EqualScoringCount = 1;

			read->clearScores(score_max_loc);
			read->Alignments = new Align[1];

			SendSeToBuffer(read, 0);
		} else {
			std::sort(read->Scores, read->Scores + read->numScores(), sortLocationScore);

			//Log.Message("%s", read->name);
//			read->EqualScoringID = 0;
			//int esid = 1;
			int n = std::min(read->numScores(), topn);
			static bool const equalOnly = false;
			int j = 1;

			if(equalOnly) {
				int i = 1;
				while(i < n && read->Scores[0].Score.f == read->Scores[i].Score.f) {
					i += 1;
				}
				n = i;
				read->EqualScoringCount = n;
			} else {
				read->EqualScoringCount = 1;
			}

			int mq = MAX_MQ;
			if(read->numScores() > 1) {
				mq = ceil(MAX_MQ * (read->Scores[0].Score.f - read->Scores[1].Score.f) / read->Scores[0].Score.f);
			}
			read->mappingQlty = mq;

			read->Alignments = new Align[n];

			for (; j < n; ++j) {
//				Log.Message("%f == %f", read->Scores[0].Score.f, read->Scores[j].Score.f);
				if(equalOnly && read->Scores[0].Score.f != read->Scores[j].Score.f) {
					Log.Error("not equal");
					Fatal();
				}
				SendSeToBuffer(read, j);
			}
			SendSeToBuffer(read, 0);
		}
	}
	else {
		read->mappingQlty = 0;
		SendSeToBuffer(read, -1);
	}

//	if (NGM.Paired() && (!read->HasFlag(NGMNames::PairedFail))) {
//		if (read->Paired == 0) {
//			Log.Error("No read pair found.");
//			Fatal();
//		}
//		Log.Verbose("[PAIRED] SW::SendToPostprocessing: %i (%s)", read->ReadId, read->name);
//		// Paired read finished score calculation
//		// -> race condition here: both paired reads could enter at the same time
//		//Log.Message("%d %d", read->Paired->Calculated, read->Paired->nScores());
//		if (read->Paired->Calculated > -1 && read->Paired->Calculated == read->Paired->numScores()) {
//			if (read->ReadId < read->Paired->ReadId) {
//				if (AtomicInc(&read->Lock) == 1) {
//					PairedReadSelection(read, read->Paired);
//				}
//			} else {
//				if (AtomicInc(&read->Paired->Lock) == 1) {
//					PairedReadSelection(read->Paired, read);
//				}
//			}
//		}
//	} else {

//}
}

int scount = 0;

void ScoreBuffer::addRead(LocationScore * newScores, int count) {

	if (count == 0) {
		Log.Error("count == 0");
		Fatal();
	}
//	//TODO: remove
//	char const * debugRead = "adb-100bp-20mio-paired.000000558.2";
//	if (strcmp(newScores[0].Read->name, debugRead) == 0) {
//		Log.Error("Adding score computations for read %s: %d", newScores[0].Read->name, count);
//		getchar();
//	}

//	if(count >= swBatchSize) {
//		Log.Warning("Read found (%d)!!", ++scount);
//	}
	for (int i = 0; i < count; ++i) {
		Log.Verbose("Adding score %d to buffer.", iScores);
		scores[iScores++] = &newScores[i];
		if(iScores == swBatchSize) {
			DoRun();
			iScores = 0;
		}
	}
}

AlignmentBuffer * out;

void ScoreBuffer::flush() {
	DoRun();
	iScores = 0;
//	NGM.bSWO.Release();
//	NGM.FinishStage(2);
	out->flush();

}

void ScoreBuffer::SendSeToBuffer(MappedRead* read, int const scoreID) {

	if (!read->hasCandidates()) {
		//NGM.AddUnmappedRead(read, MFAIL_NOCAND);
		//NGM.SaveRead(read, false);
		out->addRead(read, -1);

		Log.Message("Read %s (%i) not mapped (no candidates/scores)", read->name, read->ReadId);
		NGM.GetReadProvider()->DisposeRead(read);
		return;
	}
	//		static int const eslimit = Config.GetInt("max_equal");
	//		if (eslimit > 0 && score_max_count > eslimit) {
	//			Log.Verbose(
	//					"Read %i rejected for hitting equal scoring limit (%i top scores of %i)",
	//					read->ReadId, score_max_count, score_max);
	//
	//			score_max_count = 1;
	//			read->EqualScoringCount = 1;
	//		}

	//Log.Green("Single: %f", read->TLS()->Score.f);

	out->addRead(read, scoreID);
	//TODO: AddRead to Output
//	NGM.bSWO.Write(&read, 1);
	//NGM.AddMappedRead(read->ReadId);

//TODO: fix
//	int score_max_loc = read->TopScore;
//	float score_max = read->TLS()->Score.f;
//	int score_max_count
//	if (score_max_count > 1) {
//		Log.Error("Printing equal scoring reads is not supported at the moment!");
//		Fatal();
//		read->EqualScoringCount = score_max_count;
//
//		int esid = 1;
//		for (int j = score_max_loc + 1; j < read->nScores(); ++j) {
//			if (read->Scores[j]->Score.f == score_max) {
//
//				MappedRead * ntopread = new MappedRead(read->ReadId);
//				ntopread->EqualScoringID = esid++;
//				ntopread->EqualScoringCount = score_max_count;
//				ntopread->Calculated = 1;
//				ntopread->TopScore = 0;
//
//				ntopread->Seq = new char[qryMaxLen];
//
//				//TODO: only copy pointer, faster
//				memcpy(ntopread->Seq, read->Seq, qryMaxLen);
//				if (read->RevSeq != 0) {
//					ntopread->RevSeq = new char[qryMaxLen];
//					memcpy(ntopread->RevSeq, read->RevSeq, qryMaxLen);
//				}
//				ntopread->qlty = 0;
//				if (read->qlty != 0) {
//					ntopread->qlty = new char[qryMaxLen];
//					memcpy(ntopread->qlty, read->qlty, qryMaxLen);
//				}
//
//				ntopread->name = new char[read->nameLength];
//				memcpy(ntopread->name, read->name, read->nameLength);
//				ntopread->nameLength = read->nameLength;
//
//				ntopread->length = read->length;
//				ntopread->AddScore(*read->Scores[j])->Read = ntopread;
//				NGM.bSWO.Write(&ntopread, 1);
//
//			}
//		}
//	}

}

//void SWwoBuffer::PairedReadSelection(MappedRead * read1, MappedRead * read2) {
//	Log.Verbose("Paired Read selection -> read %i: %i candidates / read %i: %i candidates", read1->ReadId, read1->numScores(), read2->ReadId, read2->numScores());
//
//	if (read1->HasFlag(NGMNames::PairedFail) || read2->HasFlag(NGMNames::PairedFail))
//		Log.Error("Paired read selection recursion (read %i, %i)", read1->ReadId, read2->ReadId);
//
//#ifdef VERBOSE
////	for (int i = 0; i < read1->nScores(); ++i) {
////		Log.Verbose("Read %i c%i -> Location <%i, %i>, Score %f", read1->ReadId, i,
////				read1->Scores[i]->Location.m_Location, read1->Scores[i]->Location.m_RefId,
////				read1->Scores[i]->Score.f);
////	}
////
////	for (int i = 0; i < read2->nScores(); ++i) {
////		Log.Verbose("Read %i c%i -> Location <%i, %i>, Score %f", read2->ReadId, i,
////				read2->Scores[i]->Location.m_Location, read2->Scores[i]->Location.m_RefId,
////				read2->Scores[i]->Score.f); // */
////	}
//#endif
//
//	float topScore = 0;
//	int distance = 0;
//	bool equalScoreFound = false;
//
//	//Brute force
//	//TODO: fast selection for equal scoring positions.
//	int TopScore1 = -1;
//	int TopScore2 = -1;
//
//	//float maxScore1 = 0;
//	//float maxScore2 = 0;
//
//	for (int i = 0; i < read1->numScores(); ++i) {
//		//maxScore1 = std::max(maxScore1, read1->Scores[i].Score.f);
//		for (int j = 0; j < read2->numScores(); ++j) {
//			//maxScore2 = std::max(maxScore2, read2->Scores[j].Score.f);
//			if (CheckPairs(&read1->Scores[i], &read2->Scores[j], topScore, distance, equalScoreFound)) {
//				TopScore1 = i;
//				TopScore2 = j;
//				Log.Green("Accepted - %f + %f", read1->Scores[i].Score.f, read2->Scores[j].Score.f);
//			} else {
//				Log.Error("Not accepted - %f + %f", read1->Scores[i].Score.f, read2->Scores[j].Score.f);
//			}
//		}
//	}
//
//	if (TopScore1 != -1) { //&& (read1->Scores[TopScore1].Score.f + read2->Scores[TopScore2].Score.f) > (maxScore1 + maxScore2) * 0.9f) {
//		read1->TopScore = TopScore1;
//		read2->TopScore = TopScore2;
//
//		//Log.Green("Read 1: %f, Read 2: %f", read1->TLS()->Score.f, read2->TLS()->Score.f);
//
//		Log.Verbose("Read pairing found: R %i (#%i) at <%i, %i> Score %f, R %i (#%i) at <%i, %i> Score %f / Distance %i",
//				read1->ReadId, read1->TopScore, read1->Scores[read1->TopScore].Location.m_Location,
//				read1->Scores[read1->TopScore].Location.m_RefId, read1->Scores[read1->TopScore].Score.f,
//				read2->ReadId, read2->TopScore, read2->Scores[read2->TopScore].Location.m_Location,
//				read2->Scores[read2->TopScore].Location.m_RefId, read2->Scores[read2->TopScore].Score.f,
//				distance);
//
//		NGM.bSWO.Write(&read1, 1);
//		NGM.bSWO.Write(&read2, 1);
//
//		//NGM.AddMappedRead(read1->ReadId);
//		//NGM.AddMappedRead(read2->ReadId);
//
//		tCount += 1;
//		tSum += distance;
//
//		read1->clearScores();
//		read2->clearScores();
////		Log.Message("%s/%s: %d (%d)", read1->name, read2->name, distance, tSum / tCount);
//
//	} else {
//		read1->SetFlag(NGMNames::PairedFail);
//		read2->SetFlag(NGMNames::PairedFail);
//
//		Log.Verbose("Broken pair: %i, %i (%s)", read1->ReadId, read2->ReadId, read1->name);
//		brokenPairs += 1;
//
//		SendSeToBuffer(read1);
//		SendSeToBuffer(read2);
//	}
//}
//
//bool SWwoBuffer::CheckPairs(LocationScore * ls1, LocationScore * ls2, float & topScore, int & dst, bool & equalScore) {
//	if (ls1->Location.m_RefId == (ls2->Location.m_RefId ^ 0x1)) {
//		int distance =
//				(ls2->Location.m_Location > ls1->Location.m_Location) ?
//						ls2->Location.m_Location - ls1->Location.m_Location + ls2->Read->length :
//						ls1->Location.m_Location - ls2->Location.m_Location + ls1->Read->length;
//
//		if (true && distance > _NGM::sPairMinDistance && distance < _NGM::sPairMaxDistance) {
////			Log.Green("[%d, %d] Score: %f, Distance: %d (%u - %u), AVG: %d", ls1->Read->ReadId, ls2->Read->ReadId, ls1->Score.f + ls2->Score.f, distance, ls2->Location.m_Location, ls1->Location.m_Location, tSum / tCount);
//			float f = ls1->Score.f + ls2->Score.f;
//			if (f > topScore * 1.00f) {
//				topScore = f;
//				dst = distance;
//				return true;
//			} else if (f == topScore * 1.00f) {
//				equalScore = true;
//				Log.Verbose("Found a pair (%d, %d) with same score (%f). Choosing pair according to insert size.", ls1->Read->ReadId, ls2->Read->ReadId, f);
//				int avg = tSum / tCount;
//				if (abs(dst - avg) > abs(distance - avg)) {
//					topScore = f;
//					dst = distance;
//					return true;
//				} else if (abs(dst) == abs(distance)) {
//					Log.Verbose("Same insert size (%d). Choosing randomly.", distance);
//				}
//			}
//		}
//	}
//	return false;
//}
