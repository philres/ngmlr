#ifndef __ALIGNMENTDISPATCHER_H__
#define __ALIGNMENTDISPATCHER_H__

#include "IAlignment.h"
#include "Types.h"
#include "NGMThreads.h"
#include <vector>

class AlignmentDispatcher : public IAlignment
{
public:
	static AlignmentDispatcher * Instance();

	// Für den Dispatcher sind die Angaben nicht wirklich sinnvoll
	int GetScoreBatchSize() const;
	int GetAlignBatchSize() const;

	int BatchScore(
		int const mode,
		int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList,
		float * const results,
		void * extData);

	int BatchAlign(
		int const mode,
		int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList,
		Align * const results,
		void * extData);

	float GPULoad() const;
	void ResetLoad();

	void Shutdown();
private:
	static void CreateInstance();

	AlignmentDispatcher();
	~AlignmentDispatcher();

	int InitAligners();

	static AlignmentDispatcher * sInstance;

	class Aligner;

	Aligner * PickFreeAligner(int desire);
	bool ChooseCPU(int desire);
	IAlignment * GetCPUAligner();
	void FreeCPUAligner(IAlignment * aligner);

	static const int cMaxAligner = 32;

	Aligner * m_Aligner[cMaxAligner];

	NGMMutex m_SelectionMutex;
	NGMThreadWait m_SelectionWait;

	bool m_CPUEnabled;
	pfCreateAlignment m_CreateCPUAlignment;
	pfDeleteAlignment m_DeleteCPUAlignment;
	std::vector<IAlignment*> m_CPUKernels;
	NGMMutex m_CPUSelectionMutex;

	ulong m_StartTime;

	int m_TotalAligner;
	volatile int m_FreeAligner;

	int m_ScoreBatchSize;
	int m_AlignBatchSize;
};

#endif
