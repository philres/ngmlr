#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>

#include "Types.h"

#include "Log.h"
#include "Config.h"
#include "PlatformSpecifics.h"
#include "SequenceProvider.h"
#include "NGMThreads.h"
#include "MappedRead.h"


#include "NGMStats.h"

#include "NGMTask.h"

#include "IRefProvider.h"
#include "IReadProvider.h"
#include "IAlignment.h"
#include "FileWriter.h"

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

class GenericReadWriter;

#include "MapFails.h"

class _NGM
{
public:
	static _NGM & Instance();

	// Liefert den nächsten Batch an Reads, oder einen leeren Vektor wenn keien Reads mehr zu mappen sind
	std::vector<MappedRead*> GetNextReadBatch(int batchSize);
	inline bool Running() const { return m_ActiveThreads > 0; }
	inline bool DualStrand() const { return m_DualStrand; }
	inline int GetOutputFormat() const { return m_OutputFormat; }


	void GeneratePartitions();

	_NGM();
	~_NGM();

	void InitProviders();

	IAlignment * CreateAlignment(int const mode);
	void DeleteAlignment(IAlignment* instance);

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

	void * getWriter();
	void ReleaseWriter();

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

	IReadProvider * GetReadProvider();

	void FinishStage( int tid );

	NGMStats * Stats;

	static char const * AppName;
	static int sPairMinDistance;
	static int sPairMaxDistance;
private:
	static NGMTHREADFUNC ThreadFunc(void*);

	static void Init();

	friend void NGMTask::FinishStage();

	void FinishThread( int tid );
	int GetStart();
	int GetCount();

	void StartCS(int threads);

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

	//TODO: hack - fix this!!!
	void * m_Output;

	volatile int m_CurStart;
	volatile int m_CurCount;

	GenericReadWriter * writer;

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

	volatile int m_UnmappedReads;
	volatile int m_MappedReads;
	volatile int m_WrittenReads;
	volatile int m_ReadReads;


	IRefProvider * m_RefProvider;
	IReadProvider * m_ReadProvider;

	friend class _SequenceProvider;
};

#undef module_name
#define module_name 0

#define NGM _NGM::Instance()

void Help();

#endif
