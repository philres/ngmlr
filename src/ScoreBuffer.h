/*
 * SWwoBuffer.h
 *
 *  Created on: Feb 19, 2013
 *      Author: philipp_
 */

#ifndef SWWOBUFFER_H_
#define SWWOBUFFER_H_

#include "IAlignment.h"
#include "NGM.h"
#include "AlignmentBuffer.h"

#include <list>

#undef module_name
#define module_name "FILTER"

class ScoreBuffer {
public:
private:

	struct Score {
		MappedRead * read;
		int scoreId;
	};

	const int m_AlignMode;

	long pairDistCount;
	long pairDistSum;
	long brokenPairs;

private:

	bool CheckPairs(LocationScore * ls1, int const readLength1, LocationScore * ls2, int const readLength2, float & topScore, int & dst, int & equalScore);

	void top1SE(bool const equalOnly, MappedRead* read);
	void topNSE(int const topn, bool const equalOnly, MappedRead* read);
	void top1PE(bool const equalOnly, MappedRead* read);
	void topNPE(int const topn, bool const equalOnly, MappedRead* read);

	void DoRun();
	void computeMQ(MappedRead* read);

	const char** m_QryBuffer;
	const char** m_RefBuffer;
	float* m_ScoreBuffer;

	int qryMaxLen;
	int refMaxLen;

	int corridor;
	char * m_DirBuffer;
	bool m_EnableBS;

	Score * scores;
	int iScores;

	IAlignment * aligner;

	AlignmentBuffer * out;
	const int swBatchSize;

	float scoreTime;

public:
	static ulong scoreCount;

	ScoreBuffer(IAlignment * mAligner, AlignmentBuffer * mOut) :
			m_AlignMode(Config.GetInt("mode", 0, 1)), pairDistCount(1), pairDistSum(0), brokenPairs(0), aligner(mAligner), out(mOut), swBatchSize(aligner->GetScoreBatchSize() / 2) {

		m_QryBuffer = 0;
		m_RefBuffer = 0;
		m_ScoreBuffer = 0;

		corridor = Config.GetInt("corridor");

		m_QryBuffer = new char const *[swBatchSize];
		m_RefBuffer = new char const *[swBatchSize];
		m_ScoreBuffer = new float[swBatchSize];

		m_DirBuffer = new char[swBatchSize];

		m_EnableBS = false;
		//	if (Config.Exists("bs_mapping"))
		m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

		qryMaxLen = Config.GetInt("qry_max_len");
		refMaxLen = ((qryMaxLen + Config.GetInt("corridor")) | 1) + 1;

		for (int i = 0; i < swBatchSize; ++i) {
			m_RefBuffer[i] = new char[refMaxLen];
		}

		scores = new Score[swBatchSize];
		iScores = 0;
		scoreTime = 0.0f;
	}

	~ScoreBuffer() {
		Log.Verbose("SW dtor");
		delete[] m_DirBuffer;
		m_DirBuffer = 0;
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

	void flush();

	float getTime() {
		float tmp = scoreTime;
		scoreTime = 0.0f;
		return tmp;
	}

	inline int GetStage() const {
		return 2;
	}

	inline const char* GetName() const {return "SW";}
};

#endif /* SWWOBUFFER_H_ */
