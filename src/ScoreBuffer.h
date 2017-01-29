/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef SWWOBUFFER_H_
#define SWWOBUFFER_H_

#include <list>

#include "IAlignment.h"
#include "NGM.h"
#include "AlignmentBuffer.h"

#undef module_name
#define module_name "FILTER"

class ScoreBuffer {
public:
private:

	struct Score {
		MappedRead * read;
		int scoreId;
	};

private:

	void topNSE(MappedRead* read);

	void DoRun();
	void computeMQ(MappedRead* read);
	int computeMQ(float bestScore, float secondBestScore);

	void debugScoresFinished(MappedRead * read);
	ReadGroup* updateGroupInfo(MappedRead* cur_read);

	const char** m_QryBuffer;
	const char** m_RefBuffer;
	float* m_ScoreBuffer;

	int qryMaxLen;
	uloc refMaxLen;

	int corridor;

	Score * scores;
	int iScores;

	IAlignment * aligner;

	AlignmentBuffer * out;
	const int swBatchSize;

	float scoreTime;

public:

	ScoreBuffer(IAlignment * mAligner, AlignmentBuffer * mOut) :
			aligner(mAligner), out(mOut), swBatchSize(aligner->GetScoreBatchSize()) {

		m_QryBuffer = 0;
		m_RefBuffer = 0;
		m_ScoreBuffer = 0;

		corridor = Config.getReadPartCorridor();

		m_QryBuffer = new char const *[swBatchSize];
		m_RefBuffer = new char const *[swBatchSize];
		m_ScoreBuffer = new float[swBatchSize];

		qryMaxLen = Config.getReadPartLength() + 10;
		refMaxLen = ((qryMaxLen + corridor) | 1) + 1;

		for (int i = 0; i < swBatchSize; ++i) {
			m_RefBuffer[i] = new char[refMaxLen];
		}

		scores = new Score[swBatchSize];
		iScores = 0;
		scoreTime = 0.0f;

	}

	~ScoreBuffer() {
		delete[] scores;
		scores = 0;
		delete[] m_ScoreBuffer;
		m_ScoreBuffer = 0;

		for (int i = 0; i < swBatchSize; ++i) {
			delete[] m_RefBuffer[i];
			m_RefBuffer[i] = 0;
		}

		delete[] m_RefBuffer;
		m_RefBuffer = 0;
		delete[] m_QryBuffer;
		m_QryBuffer = 0;
	}

	void addRead(MappedRead * read, int count);

	void scoreShortRead(MappedRead * read);

	void flush();

	float getTime() {
		float tmp = scoreTime;
//		scoreTime = 0.0f;
		return tmp;
	}

	inline int GetStage() const {
		return 2;
	}

	inline const char* GetName() const {return "SW";}
};

#endif /* SWWOBUFFER_H_ */
