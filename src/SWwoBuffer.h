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
#include "Output.h"

#include <list>

#undef module_name
#define module_name "SWwoBuffer"

class SWwoBuffer {
public:
private:
	const int m_AlignMode;
//	int m_Ivk_mode;
//	int m_Ivk_batchSize;
//	const char* const * m_Ivk_refSeqList;
//	const char* const * m_Ivk_qrySeqList;
//	Align* m_Ivk_results;
//	int m_Ivk_return;
	void LaunchInvoked();

public:
	void SendToPostprocessing(MappedRead* read);
	void SendSeToBuffer(MappedRead* read);

private:
	static void PairedReadSelection(MappedRead* read1, MappedRead* read2);
	static bool CheckPairs(LocationScore* ls1, LocationScore* ls2, float&, int&, bool& equlScore);

	const char** m_QryBuffer;
	const char** m_RefBuffer;
	float* m_ScoreBuffer;

	int qryMaxLen;
	int refMaxLen;

	int corridor;
	//int batchSize;
	char * m_DirBuffer;
	bool m_EnableBS;

	LocationScore * * scores;
	int iScores;

	IAlignment * aligner;
	Output * out;
	const int swBatchSize;

	//std::list<MappedRead*> m_ReadBuffer;
	//	bool CommitReads();
public:
	static ulong scoreCount;

	SWwoBuffer(IAlignment * mAligner, Output * mOut) :
			m_AlignMode(Config.GetInt("mode", 0, 1)), aligner(mAligner), out(mOut), swBatchSize(aligner->GetScoreBatchSize() / 2) {

		m_QryBuffer = 0;
		m_RefBuffer = 0;
		m_ScoreBuffer = 0;

		corridor = Config.GetInt("corridor");

//		swBatchSize = 4096;

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

		scores = new LocationScore * [swBatchSize];
		iScores = 0;
	}

	~SWwoBuffer() {
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

	void addRead(LocationScore * scores, int count);

	void DoRun();

	void flush();

	inline int GetStage() const {
		return 2;
	}

	inline const char* GetName() const {return "SW";}
};

#endif /* SWWOBUFFER_H_ */
