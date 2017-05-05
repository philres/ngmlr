/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "ScoreBuffer.h"

#include <algorithm>
#include <cmath>

#include "NGMTask.h"
#include "Timing.h"
#include "CS.h"
#include "AlignmentBuffer.h"
#include "StrippedSW.h"

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

bool sortLocationScoreLocation(LocationScore a, LocationScore b) {
	return a.Location.m_Location < b.Location.m_Location;
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
		delete[] tmpScores;
		tmpScores = 0;
	}
}

ReadGroup* ScoreBuffer::updateGroupInfo(MappedRead* cur_read) {
	ReadGroup* group = cur_read->group;
	if (group != 0) {
		//TODO: make atomic, parts from a group can end up in different threads!
		group->readsFinished += 1;
		//				Log.Message("Scorecount for %s (%d): %d - %d", cur_read->name, cur_read->ReadId, cur_read->numScores(), cur_read->Calculated);
		if (cur_read->Scores[0].Location.isReverse()) {
			group->reverseMapped += 1;
		} else {
			group->fwdMapped += 1;
		}
		group->bestScoreSum += (int) (cur_read->Scores[0].Score.f);
	}
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
			if (loc.isReverse()) {
				// RefId auf +-Strang setzen
				//--loc.m_RefId;
				cur_read->computeReverseSeq();
				m_QryBuffer[i] = cur_read->RevSeq;
			} else {
				m_QryBuffer[i] = cur_read->Seq;
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
		int res = aligner->BatchScore(0, iScores, m_RefBuffer, m_QryBuffer, m_ScoreBuffer, 0);

		Log.Debug(16, "INFO\tSCORES\t%d scores computed (out of %d)", res, iScores);

		if (res != iScores)
		Log.Error("Kernel couldn't calculate all scores (%i out of %i)", res, iScores);

		//Process results
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
				if(group != 0) {
					//If all reads from group are finished
					if(group->readsFinished == group->readNumber) {
						out->processLongReadLIS(group);
//					out->WriteRead(group->fullRead, false);
					}
				} else {
					out->processShortRead(cur_read);
				}

			}
		}
		scoreTime += tmr.ET();
	} else {
		Log.Debug(16, "INFO\tSCORES\tEmpty buffer submitted.");
	}
}

void ScoreBuffer::topNSE(MappedRead* read) {

	//Sort scores
	std::sort(read->Scores, read->Scores + read->numScores(),
			sortLocationScore);

	int numScores = read->numScores();

	if(numScores > 1) {
		float const minScore = read->Scores[0].Score.f * 0.75f;
		int i = 1;
		while(i < numScores && read->Scores[i].Score.f > minScore) {
			i += 1;
		}
		numScores = i;
	}

	//compute mapping quality
	computeMQ(read);

	//numScores alignments will be computed in the next step
	read->Calculated = numScores;
}


void ScoreBuffer::addRead(MappedRead * read, int count) {

	LocationScore * newScores = read->Scores;
	if (count == 0) {
		Log.Error("Internal error (count == 0). Please report this on https://github.com/Cibiv/NextGenMap/issues");
	}

	//Adding scores to buffer. If buffer full, submit to CPU/GPU for score computation
	for (int i = 0; i < count; ++i) {

		Log.Debug(256, "READ_%d\tSCORES_BUFFER\tCMR_%d %f (location %llu) added to score buffer at position %d", read->ReadId, i, newScores[i].Score.f, newScores[i].Location.m_Location, iScores);

		scores[iScores].read = read;
		scores[iScores++].scoreId = i;
		if(iScores == swBatchSize) {
			DoRun();
			iScores = 0;
		}
	}
}

void ScoreBuffer::scoreShortRead(MappedRead * read) {
	IAlignment * aligner = new StrippedSW();

	static int const readPartLength = Config.getReadPartLength();

	// Remove redundant candidates (CMRs in close proximity)
	LocationScore * tmpScore = new LocationScore[read->numScores()];
	int tmpScoreIndex = 0;

	std::sort(read->Scores, read->Scores + read->numScores(), sortLocationScoreLocation);

	uloc lastLocation = 0;

	for (int i = 0; i < read->numScores(); ++i) {
		if(lastLocation - read->Scores[i].Location.m_Location > readPartLength) {
			tmpScore[tmpScoreIndex++] = read->Scores[i];
		}
		lastLocation = read->Scores[i].Location.m_Location;
	}

	delete[] read->Scores;
	read->Scores = 0;

	read->AllocScores(tmpScore, tmpScoreIndex);

	delete[] tmpScore;
	tmpScore = 0;


	for (int i = 0; i < read->numScores(); ++i) {
		int corridor = read->length * 0.3 + 256;

		char * qrySeq = 0;
		char * refSeq = new char[read->length + corridor + 10];
		memset(refSeq, '\0', read->length + corridor + 10);
		if (read->Scores[i].Location.isReverse()) {
			read->computeReverseSeq();
			qrySeq = read->RevSeq;
		} else {
			qrySeq = read->Seq;
		}

		//decode reference sequence
		if (!SequenceProvider.DecodeRefSequence(refSeq, 0,
				read->Scores[i].Location.m_Location - (corridor >> 1), read->length + corridor)) {
			//							loc.m_Location - (corridor >> 1), cur_read->length + corridor)) {
//			Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", loc.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
			//Log.Warning("Read sequence: %s", cur_read->Seq);
			memset(refSeq, 'N', read->length + corridor);
		}
		float score = 0.0f;
		aligner->SingleScore(0, corridor, refSeq, qrySeq, score, 0);
		read->Scores[i].Score.f = score;

		delete[] refSeq;
		refSeq = 0;
	}


	std::sort(read->Scores, read->Scores + read->numScores(), sortLocationScore);

//	Log.Message("Read: %s", read->name);
//	for (int i = 0; i < read->numScores(); ++i) {
//		Log.Message("Score %d: %f", i, read->Scores[i].Score.f);
//	}
	computeMQ(read);

	out->processShortRead(read);

	delete aligner;
}

void ScoreBuffer::flush() {
	//Force submitting remaining computation from buffer to CPU/GPU
	DoRun();
	iScores = 0;
//	out->flush();
}

