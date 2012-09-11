#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>

#include "Types.h"

#include "Log.h"
#include "Config.h"
#include "PlatformSpecifics.h"
#include "SequenceProvider.h"
#include "Buffer.h"
#include "NGMThreads.h"
#include "MappedRead.h"

#include "NGMStats.h"

#include "NGMTask.h"

#include "Aligner.h"
//#include "Partition.h"
#include "IRefProvider.h"
#include "IReadProvider.h"

//#include "ReadBuffer.h"
//#include "CSCache.h"

#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifndef NDEBUG
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
#endif

#undef module_name
#define module_name "NGM"

#include "MapFails.h"

class _NGM
{
public:
	static _NGM & Instance();

	// Liefert den nächsten Batch an Reads, oder einen leeren Vektor wenn keien Reads mehr zu mappen sind
	std::vector<MappedRead*> GetNextReadBatch(int batchSize);
	inline bool Running() const { return m_ActiveThreads > 0; }
	inline bool DualStrand() const { return m_DualStrand; }
	inline bool Paired() const { return m_Paired; }
	inline int GetOutputFormat() const { return m_OutputFormat; }

	IAlignment * Aligner() const { return AlignmentDispatcher::Instance(); }

	void GeneratePartitions();

	_NGM();
	~_NGM();

	void InitProviders();

	void StartThread( NGMTask * task, int cpu = -1 );

	void AddUnmappedRead(MappedRead const * const read, int reason );
	int GetUnmappedReadCount() const;
	void AddMappedRead( int readid );
	int GetMappedReadCount() const;
	void AddWrittenRead( int readid );
	int GetWrittenReadCount() const;
	void AddReadRead( int readid );
	int GetReadReadCount() const;

	void DisposeRead(MappedRead * read);

	void AquireOutputLock();
	void ReleaseOutputLock();

	void InitQuit();
	void StartThreads();
	void MainLoop();

	bool StageActive( int stage )
	{
		return GetStageThreadCount(stage) > 0;
	}
	bool SingleStageActive( int stage )
	{
		return m_StageThreadCount[stage] > 0;
	}
	int GetStageThreadCount( int stage );

	bool ThreadActive( int tid, int stage );

	IRefProvider const * GetRefProvider(int const tid);
	void ReleaseRefProvider(int const tid);

	IReadProvider * GetReadProvider();
	void ReleaseReadProvider();

	void SaveRead(MappedRead* read, bool mapped = true);


	//bool LastPartition();
	//void BufferRead(MappedRead * read);

	NGMStats * Stats;

	// Switch to a list accessed by Stage?
	Buffer<LocationScore*> bCSSW;	// Buffer CS -> SW
	Buffer<MappedRead*> bSWO;	// Buffer SW -> Output

	static char const * AppName;
	static int sPairMinDistance;
	static int sPairMaxDistance;
private:
	static NGMTHREADFUNC ThreadFunc(void*);

	static void Init();

	friend void NGMTask::FinishStage();
	void FinishStage( int tid );
	void FinishThread( int tid );
	int GetStart();
	int GetCount();

	void StartCS(int threads);
	void StartSW(int threads);

	void UpdateScheduler(float, float);
	bool CanSwitch();

	friend int main(int argc, char* argv[]);
	static _NGM * pInstance;
	static NGMOnceControl once_control;
	static const int cMaxThreads = 1024;
	static const int cMaxStage = 8;

	volatile int m_ActiveThreads;
	int m_NextThread;
	bool const m_DualStrand;
	bool const m_Paired;
	int const m_OutputFormat;
	int m_ReadStart;
	int m_ReadCount;

	volatile int m_CurStart;
	volatile int m_CurCount;
//	volatile int m_CurrentPartition;
//	volatile int m_ReadsBuffered;

//	volatile bool m_RefGenPending;

	NGMMutex m_Mutex;
	NGMMutex m_OutputMutex;
	int m_StageThreadCount[cMaxStage];
	NGMThread m_Threads[cMaxThreads];
	NGMTask * m_Tasks[cMaxThreads];
	int m_BlockedThreads[cMaxStage];
	int m_ToBlock[cMaxStage];
	NGMMutex m_SchedulerMutex;
	NGMThreadWait m_SchedulerWait;
	NGMMutex m_UMRMutex;
	NGMThreadWait m_CSWait;
	bool m_TrackUnmappedReads;
//	std::vector<int> m_UnmappedReadList;
	volatile int m_UnmappedReads;
	volatile int m_MappedReads;
	volatile int m_WrittenReads;
	volatile int m_ReadReads;

	//PartitionList m_Partitions;
	//IRefProvider * m_CurrentRef;
	IRefProvider * m_RefProvider;
	IReadProvider * m_ReadProvider;
	//CSCache * m_Cache;

	//ReadBuffer * m_ReadBuffer;

//	PartitionList GeneratePartitions(std::list<int> SeqIds, long PartitionThreshhold, long Overlap);
//	std::list<int> GetRefConfig();

	void UpdateReadCount(int n);

	friend class _SequenceProvider;
};

#undef module_name
#define module_name 0

#define NGM _NGM::Instance()

void Fatal();
void Help();

#endif
