#include "NGM.h"

#include <memory.h>
#include <stdlib.h>

#include "Debug.h"
#include "PrefixTable.h"
#include "ReadProvider.h"
#include "FileWriter.h"
#include "SAMWriter.h"
#include "BAMWriter.h"

#include <limits.h>

#include "CS.h"

#undef module_name
#define module_name "NGM"

_NGM * _NGM::pInstance = 0;
NGMOnceControl _NGM::once_control = NGM_ONCE_INIT;
char const * _NGM::AppName = 0;
int _NGM::sPairMinDistance = 0;
int _NGM::sPairMaxDistance = INT_MAX;

namespace __NGM {
inline int min(int a, int b) {
	return (a < b) ? a : b;
}
}

void _NGM::Init() {
	pInstance = new _NGM();

//	sPairMinDistance = Config.GetInt("pair_min_distance", 0, -1);
//	sPairMaxDistance = Config.GetInt("pair_max_distance", -1, -2);
	sPairMinDistance = Config.GetInt("min_insert_size");
	sPairMaxDistance = Config.GetInt("max_insert_size");
	if (sPairMaxDistance <= 0)
		sPairMaxDistance = INT_MAX;

}

_NGM & _NGM::Instance() {
	NGMOnce(&_NGM::once_control, Init);
	return *pInstance;
}

_NGM::_NGM() :
		Stats(NGMStats::InitStats(AppName)), m_ActiveThreads(0), m_NextThread(0), m_DualStrand(Config.GetInt("dualstrand") != 0), m_Paired(
		Config.GetInt("paired") != 0 || (Config.Exists("qry1") && Config.Exists("qry2"))),
#ifdef _BAM
				m_OutputFormat(Config.GetInt("format", 0, 2)),
#else
				m_OutputFormat(Config.GetInt("format", 0, 1)),
#endif
				m_CurStart(0), m_CurCount(0), m_SchedulerMutex(), m_SchedulerWait(), m_TrackUnmappedReads(false), m_UnmappedReads(0), m_MappedReads(
						0), m_WrittenReads(0), m_ReadReads(0), m_ReadProvider(0) {

	char const * const output_name = Config.GetString("output");
	if (m_OutputFormat != 2) {
		Log.Message("Opening for output (SAM): %s", output_name);
		m_Output = new FileWriter(output_name);
	} else {
		Log.Message("Opening for output (BAM): %s", output_name);
		m_Output = new FileWriterBam(output_name);
	}

	Log.Message("NGM Core initialization");
	NGMInitMutex(&m_Mutex);
	NGMInitMutex(&m_OutputMutex);
	NGMInitMutex(&m_UMRMutex);
	NGMInitWait(&m_CSWait);
	NGMInitMutex(&m_SchedulerMutex);
	NGMInitWait(&m_SchedulerWait);

	writer = 0;

	memset(m_StageThreadCount, 0, cMaxStage * sizeof(int));
	memset(m_BlockedThreads, 0, cMaxStage * sizeof(int));
	memset(m_ToBlock, 0, cMaxStage * sizeof(int));
	if (m_Paired && !m_DualStrand)
		Log.Error("Logical error: Paired read mode without dualstrand search.");

	}

void _NGM::InitProviders() {
	CS::Init();

	SequenceProvider.Init(); // Prepares input data

	m_RefProvider = new CompactPrefixTable();

	if (Config.Exists("qry") || (Config.Exists("qry1") && Config.Exists("qry2"))) {
		m_ReadProvider = new ReadProvider();
		uint readCount = m_ReadProvider->init();
//		Log.Message("Read count %d", readCount);
//		NGM.UpdateReadCount(readCount);
	}
}

_NGM::~_NGM() {
	//if (m_ReadBuffer != 0)
	//	delete m_ReadBuffer;

	if (m_RefProvider != 0)
		delete m_RefProvider;

	if (m_ReadProvider != 0)
		delete m_ReadProvider;
}

int _NGM::GetStart() {
	int start = Config.GetInt("qry_start");
	if (start < 0) {
		Log.Error("Invalid start read %i", start);
		start = 0;
	} else {
		start *= ((m_Paired) ? 2 : 1);
	}
	return start;
}
int _NGM::GetCount() {
	int count = Config.GetInt("qry_count");
	if (count > 0)
		count *= ((m_Paired) ? 2 : 1);
	return count;
}

void _NGM::StartThread(NGMTask * task, int cpu) {
	Log.Verbose("Starting thread %i <%s> on cpu %i", m_NextThread, task->GetName(), cpu);
	NGMLock(&m_Mutex);
	task->m_TID = m_NextThread;
	m_Tasks[m_NextThread] = task;
	m_Threads[m_NextThread] = NGMCreateThread(&_NGM::ThreadFunc, task, true);

	if (cpu != -1)
	NGMSetThreadAffinity(&m_Threads[m_NextThread], cpu);

	++m_StageThreadCount[task->GetStage()];
	++m_NextThread;
	++m_ActiveThreads;
	NGMUnlock(&m_Mutex);
}

void _NGM::AddUnmappedRead(MappedRead const * const read, int reason) {
	AtomicInc(&m_UnmappedReads);
//	m_Output->SaveRead(read, false);
	Log.Verbose("Read %s (%i) not mapped (%i)", read->name, read->ReadId, reason);
//	if (m_TrackUnmappedReads) {
//		NGMLock(&m_UMRMutex);
//		m_UnmappedReadList.push_back(readid);
//		NGMUnlock(&m_UMRMutex);
//	}
}

void * _NGM::getWriter() {
	return m_Output;
}

void _NGM::ReleaseWriter() {
	if (m_Output != 0) {
		if (m_OutputFormat != 2) {
			delete (FileWriter*) m_Output;
		} else {
			delete (FileWriterBam*) m_Output;
		}
		m_Output = 0;
	}
}

int _NGM::GetUnmappedReadCount() const {
	return m_UnmappedReads;
}

void _NGM::AddMappedRead(int readid) {
	AtomicInc(&m_MappedReads);
}
int _NGM::GetMappedReadCount() const {
	return m_MappedReads;
}

void _NGM::AddWrittenRead(int readid) {
	AtomicInc(&m_WrittenReads);
}
int _NGM::GetWrittenReadCount() const {
	return m_WrittenReads;
}

void _NGM::AddReadRead(int readid) {
	AtomicInc(&m_ReadReads);
}

int _NGM::GetReadReadCount() const {
	return m_ReadReads;
}

int _NGM::GetStageThreadCount(int stage) {
	NGMLock(&m_Mutex);
	int cnt = 0;
	for (int i = 0; i <= stage; ++i) {
		cnt += m_StageThreadCount[i];
	}
	NGMUnlock(&m_Mutex);
	return cnt;
}

void _NGM::FinishStage(int tid) {
	NGMLock(&m_Mutex);
	int stage = m_Tasks[tid]->GetStage();
	--m_StageThreadCount[stage];
	NGMUnlock(&m_Mutex);
	Log.Verbose("Thread %i finished its stage (Stage %i, now %i active)", tid, stage, m_StageThreadCount[stage]);
}

void _NGM::FinishThread(int tid) {
	AtomicDec(&m_ActiveThreads);

	Log.Verbose("Thread %i finished (%i worker threads remaining)", tid, m_ActiveThreads);
	m_Tasks[tid]->FinishStage();
	delete m_Tasks[tid];
	m_Tasks[tid] = 0;
}

//void _NGM::UpdateReadCount(int n) {
//	if ((m_ReadCount < 0) || ((m_ReadStart + m_ReadCount) > n))
//		m_ReadCount = n - m_ReadStart;
//
//	/*	if (Config.Exists("continue_from"))
//	 {
//	 int bufferPartition = -1;
//	 m_ReadBuffer = ReadBuffer::ContinueFrom(std::string(Config.GetString("continue_from")), m_ReadStart, m_ReadCount, bufferPartition);
//	 m_CurrentPartition = bufferPartition;
//	 m_ReadsBuffered = m_ReadCount;
//	 Log.Green("Restarting from partition %i", m_CurrentPartition);
//	 }
//	 else
//	 {*/
//	//m_ReadBuffer = ReadBuffer::GenerateNew(m_ReadStart, m_ReadCount);
//	//}*/
//}

#ifdef _DEBUG
void TestMem();
#endif

bool eof = false;

std::vector<MappedRead*> _NGM::GetNextReadBatch(int desBatchSize) {
	NGMLock(&m_Mutex);

	std::vector<MappedRead*> list;

	if (m_Paired) {
		if (desBatchSize == 1)
			desBatchSize = 2;
		else
			desBatchSize &= ~1;
	}

	if (m_CurCount == 0) {
//		m_CurStart = m_ReadStart;
//		m_CurCount = m_ReadCount;
		m_CurStart = 0;
		NGMSignal(&m_CSWait);
	}

//	if (desBatchSize > m_CurCount)
//		desBatchSize = m_CurCount;

	list.reserve(desBatchSize);
	int count = 0;

	for (int i = 0; i < desBatchSize && !eof; i = i + 2) {
		MappedRead * read1 = 0;
		MappedRead * read2 = 0;
		eof = !NGM.GetReadProvider()->GenerateRead(m_CurStart + i, read1, m_CurStart + i + 1, read2);

		if (read1 != 0) {
			count += 1;
			list.push_back(read1);
			if (read2 != 0) {
				count += 1;
				list.push_back(read2);
			}
		}
	}

	m_CurStart += count;
	//m_CurStart += desBatchSize;
	m_CurCount -= desBatchSize;

	NGMUnlock(&m_Mutex);

#ifdef _DEBUGCS
	if(m_CurStart > 100000) {
		Log.Warning("Debug CS mode: quitting after 100000 reads!");
		list.clear();
	}
#endif
	return list;
}

//int GetTopN() {
//	if (Config.Exists("topn"))
//		return Config.GetInt("topn", 1, -1);
//	else
//		return 10;
//}

NGMTHREADFUNC _NGM::ThreadFunc(void* data) {
	NGMTask * task = (NGMTask*) data;
	int tid = task->m_TID;

//	try {
	Log.Verbose("Running thread %i", tid);

	task->Run();

	Log.Verbose("Thread %i run return, finishing", tid);

	NGM.FinishThread(tid);

	Log.Verbose("Thread %i finished", tid);
//	} catch (...) {
//		Log.Error("Unhandled exception in thread %i", tid);
//	}

	//Log.Message("ThreadFunc on thread %i returning", tid);

	return 0;
}

void _NGM::InitQuit() {
	static int quitState = 0;
	++quitState;
	if (quitState == 1) {
		Log.Warning("Hit 'Q' two more times to quit program.");
	} else if (quitState >= 3) {
		CleanupPlatform();
		Log.Error("Terminate by user request");
		Log.Message("%i Threads still active", m_ActiveThreads);
		for (int i = 0; i < cMaxStage; ++i) {
			if (m_StageThreadCount[i] > 0)
			Log.Message("Stage %i got %i threads still running", i, m_StageThreadCount[i]);
		}
		exit(-1);
	}
}

void _NGM::AquireOutputLock() {
	NGMLock(&m_OutputMutex);
}
void _NGM::ReleaseOutputLock() {
	NGMUnlock(&m_OutputMutex);
}

IRefProvider const * _NGM::GetRefProvider(int const tid) {
	return m_RefProvider;
}

IReadProvider * _NGM::GetReadProvider() {
	return m_ReadProvider;
}

//void _NGM::ReleaseRefProvider(int const tid) {
//
//}
//
//void _NGM::ReleaseReadProvider() {
//
//}

bool _NGM::ThreadActive(int tid, int stage) {
	if (m_ToBlock[stage] > 0) {
		NGMLock(&m_SchedulerMutex);
		bool blocked = false;

		if (m_ToBlock[stage] > 0) {
			--m_ToBlock[stage];
			blocked = true;
			m_BlockedThreads[stage]++;

			while (blocked) {
				Log.Green("Block %i @ %i", tid, stage);
				NGMWait(&m_SchedulerMutex, &m_SchedulerWait);
				if (m_ToBlock[stage] < 0) {
					++m_ToBlock[stage];
					blocked = false;
					m_BlockedThreads[stage]--;
					Log.Green("Unblock %i @ %i", tid, stage);
				}
			}
		}
		NGMUnlock(&m_SchedulerMutex);
	}
	return true;
}

//// Called every 50ms
//void _NGM::UpdateScheduler(float load1, float load2) {
//	return; // disabled
//	if (load1 > 0.75f) {
//		NGMLock(&m_SchedulerMutex);
//		if (m_BlockedThreads[2] > 1 && CanSwitch()) {
//			++m_ToBlock[0];
//			--m_ToBlock[2];
//			NGMSignal(&m_SchedulerWait);
//			Log.Green("Switching thread CS->SW");
//		}
//		NGMUnlock(&m_SchedulerMutex);
//	} else if (load1 < 0.25f) {
//		NGMLock(&m_SchedulerMutex);
//		if (m_BlockedThreads[0] > 1 && CanSwitch()) {
//			--m_ToBlock[0];
//			++m_ToBlock[2];
//			NGMSignal(&m_SchedulerWait);
//			Log.Green("Switching thread SW->CS");
//		}
//		NGMUnlock(&m_SchedulerMutex);
//	}
//}

//bool _NGM::CanSwitch() {
//	static ulong lastswitch = 0;
//	ulong now = GetTickCount();
//	if (lastswitch + 1000 < now) {
//		lastswitch = now;
//		return true;
//	} else
//		return false;
//}

void _NGM::StartThreads() {
	int threadcount = Config.GetInt("cpu_threads", 1, 0);

#ifdef _DEBUGCS
	threadcount = 1;
#endif
	//m_ToBlock[0] = 1;
	//m_ToBlock[2] = threadcount-1;

	StartCS(threadcount);
	//int swthreads = threadcount / 1;
	//int swthreads = threadcount;
	//if (swthreads < 2)
	//	swthreads = 2;
	//int swthreads = 1;
	//StartSW(1);

}

void _NGM::StartCS(int cs_threadcount) {
	int * cpu_affinities = new int[cs_threadcount];
	bool fixCPUs = Config.HasArray("cpu_threads");

	if (fixCPUs)
		Config.GetIntArray("cpu_threads", cpu_affinities, cs_threadcount);

	for (int i = 0; i < cs_threadcount; ++i) {
		NGMTask * cs = 0;
		cs = new CS();
		StartThread(cs, (fixCPUs) ? cpu_affinities[i] : -1);
	}
	delete[] cpu_affinities;
}

//void _NGM::StartSW(int sw_threadcount) {
//	for (int i = 0; i < sw_threadcount; ++i) {
//		SW * sw = new SW();
//		NGM.StartThread(sw);
//	}

//m_Output = new Output(Config.GetString("output"));
//NGM.StartThread(m_Output);
//}

#include "OclHost.h"
#include "SWOclCigar.h"
#include "seqan/EndToEndAffine.h"

IAlignment * _NGM::CreateAlignment(int const mode) {
	//int dev_type = Config.GetInt("ocl_device");
	int dev_type = CL_DEVICE_TYPE_CPU;

	if (Config.Exists("gpu")) {
		dev_type = CL_DEVICE_TYPE_GPU;
	}

	Log.Verbose("Mode: %d GPU: %d", mode, mode & 0xFF);

	OclHost * host = new OclHost(dev_type, mode & 0xFF, Config.GetInt("cpu_threads"));

	IAlignment * instance = 0;

//#ifndef NDEBUG
	//Log.Error("Alignment mode: %d", mode);
//#endif
	int ReportType = (mode >> 8) & 0xFF;
	switch (ReportType) {
		case 0:

			Log.Verbose("Output: text");

//		instance = new SWOclAlignment(host);
			//			instance = new SWOclCigar(host);
			break;
			case 1:

			Log.Verbose("Output: cigar");

			instance = new SWOclCigar(host);
			break;
			default:
			Log.Error("Unsupported report type %i", mode);
			break;
		}
	instance = new EndToEndAffine();
	return instance;
}

void _NGM::DeleteAlignment(IAlignment* instance) {
	SWOcl * test = (SWOcl *) instance;
	OclHost * host = test->getHost();
#ifndef NDEBUG
	Log.Message("Delete alignment called");
#endif
	if (instance != 0) {
		delete instance;
		instance = 0;
	}

	if (host != 0) {
		//delete host;
		host = 0;
	}
}

//TODO: remove
#include "MappedRead.h"
#include "LocationScore.h"

volatile bool Terminating = false;

void _NGM::MainLoop() {

	bool const isPaired = Config.GetInt("paired") > 0;
	float loads[2] = { 0, 0 };
	//_Log::FilterLevel(1);
#ifdef _WIN32
	char const * const loadbar = "\xFE\xFE\xFE\xFE\xFE\xFE\xFE\xFE\xFE\xFE";
#else
	//char const * const loadbar = "##########";
#endif
#ifdef INSTANCE_COUNTING
	int count = 0;
#endif
	int const threadcount = Config.GetInt("cpu_threads", 1, 0);
	bool const progress = Config.GetInt("no_progress") != 1;
	while (Running()) {
		Sleep(50);
		while (!Terminating && _kbhit()) {
			char ch = _getch();
			if (ch == 'q')
				NGM.InitQuit();
			}
			//loads[0] *= 0.25;
			//loads[0] += (NGM.bCSSW.Load() * 0.75f);
			//loads[1] /= 2;
			//loads[1] += (NGM.bSWO.Load() / 2);

//		UpdateScheduler(loads[0], loads[1]);
		if (progress) {
			int processed = std::max(1, NGM.GetMappedReadCount() + NGM.GetUnmappedReadCount());
			if (!isPaired) {
				Log.Progress("Mapped: %d, CRM/R: %d, CS: %d (%d), R/S: %d, Time: %.2f %.2f %.2f", processed, NGM.Stats->avgnCRMS, NGM.Stats->csLength, NGM.Stats->csOverflows, NGM.Stats->readsPerSecond * threadcount, NGM.Stats->csTime, NGM.Stats->scoreTime, NGM.Stats->alignTime);
				//Log.Progress("%d/%d (%.2f%), Buffer: %.2f %.2f, CS: %d %d %.2f", processed, m_ReadCount, processed * 100.0f / m_ReadCount, loads[0], loads[1], NGM.Stats->csLength, NGM.Stats->csOverflows, NGM.Stats->csTime);
			} else {
				Log.Progress("Mapped: %d, CRM/R: %d, CS: %d (%d), R/S: %d, Time: %.2f %.2f %.2f, Pairs: %.2f %.2f", processed, NGM.Stats->avgnCRMS, NGM.Stats->csLength, NGM.Stats->csOverflows, NGM.Stats->readsPerSecond * threadcount, NGM.Stats->csTime, NGM.Stats->scoreTime, NGM.Stats->alignTime, NGM.Stats->validPairs, NGM.Stats->insertSize);
				//Log.Progress("%d/%d (%.2f%), Buffer: %.2f %.2f, CS: %d %d %.2f, Pairs: %.2f %.2f", processed, m_ReadCount, processed * 100.0f / m_ReadCount, loads[0], loads[1], NGM.Stats->csLength, NGM.Stats->csOverflows, NGM.Stats->csTime, NGM.Stats->brokenPairs, NGM.Stats->insertSize);
			}
		}
#ifdef INSTANCE_COUNTING
		if (count++ == 10) {
			Log.Warning("MappedRead count = %i (max %i) (%d)", MappedRead::sInstanceCount, MappedRead::maxSeqCount, sizeof(MappedRead));
			Log.Warning("LocationScore count = %i (%d)", LocationScore::sInstanceCount, sizeof(LocationScore));
			count = 0;
		}
#endif
	}
}

