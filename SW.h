#ifndef __SW_H__
#define __SW_H__

#include "NGM.h"
#include "NGMThreads.h"

#include "IAlignment.h"

#include <list>

#undef module_name
#define module_name "SW"

class SW : public NGMTask
{
private:
	const int m_AlignMode;
	int m_Ivk_mode;
	int m_Ivk_batchSize;
	const char* const * m_Ivk_refSeqList;
	const char* const * m_Ivk_qrySeqList;
	Align* m_Ivk_results;
	int m_Ivk_return;
	void LaunchInvoked();

public:
	static void SendToPostprocessing(MappedRead* read);
	static void SendSeToBuffer(MappedRead* read);

private:
	static void PairedReadSelection(MappedRead* read1, MappedRead* read2);
	static bool CheckPairs(LocationScore* ls1, LocationScore* ls2, float&, int&, bool& equlScore);

	const char** m_QryBuffer;
	const char** m_RefBuffer;
	float* m_ScoreBuffer;

	//std::list<MappedRead*> m_ReadBuffer;
	//	bool CommitReads();
public:
	static ulong scoreCount;

	SW() :
			m_AlignMode(Config.GetInt("mode", 0, 1)) {
		NGM.bSWO.Register();
		m_QryBuffer = 0;
		m_RefBuffer = 0;
		m_ScoreBuffer = 0;
	}

	~SW() {
		Log.Verbose("SW dtor");
	}

	void DoRun();

	inline int GetStage() const {
		return 2;
	}

	inline const char* GetName() const { return "SW"; }
};

#endif
