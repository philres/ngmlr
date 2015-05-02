#include "NGM.h"

#include <memory.h>
#include <stdlib.h>
#include <limits.h>

#include "CS.h"
#include "Debug.h"
#include "PrefixTable.h"
#include "ReadProvider.h"
#include "PlainFileWriter.h"
#include "GZFileWriter.h"
#include "SAMWriter.h"
#include "BAMWriter.h"

#include "OclHost.h"
#include "SWOclCigar.h"
#include "seqan/EndToEndAffine.h"

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
		Stats(new NGMStats()), m_ActiveThreads(0), m_NextThread(0), m_DualStrand(Config.GetInt("dualstrand") != 0), m_Paired(
		Config.GetInt("paired") != 0 || (Config.Exists("qry1") && Config.Exists("qry2"))),
#ifdef _BAM
				m_OutputFormat(Config.GetInt("format", 0, 2)),
#else
				m_OutputFormat(Config.GetInt("format", 0, 1)),
#endif
				m_CurStart(0), m_CurCount(0), m_SchedulerMutex(), m_SchedulerWait(), m_TrackUnmappedReads(false), m_UnmappedReads(0), m_MappedReads(
						0), m_WrittenReads(0), m_ReadReads(0), m_ReadProvider(0) {

	char const * const output_name = Config.Exists("output") ? Config.GetString("output") : 0;
	if (m_OutputFormat != 2) {
		if (output_name != 0) {
			Log.Message("Opening for output (SAM): %s", output_name);
		} else {
			Log.Message("Wrinting output (SAM) to stdout");
		}
		if(Config.Exists(GZ)) {
			m_Output = new GZFileWriter(output_name);
		} else {
			m_Output = new PlainFileWriter(output_name);
		}
	} else {
		if (output_name != 0) {
			Log.Message("Opening for output (BAM): %s", output_name);
		} else {
			Log.Message("Wrinting output (BAM) to stdout");
		}
		m_Output = new FileWriterBam(output_name);
	}

	Log.Verbose("NGM Core initialization");
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

	m_RefProvider = new CompactPrefixTable(NGM.DualStrand());

	if (Config.Exists("qry") || (Config.Exists("qry1") && Config.Exists("qry2"))) {
		m_ReadProvider = new ReadProvider();
		uint readCount = m_ReadProvider->init();
	}
}

_NGM::~_NGM() {

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
	m_Threads[m_NextThread] = NGMCreateThread(&_NGM::ThreadFunc, task, false);

	++m_StageThreadCount[task->GetStage()];
	++m_NextThread;
	++m_ActiveThreads;
	NGMUnlock(&m_Mutex);
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

void _NGM::AddUnmappedRead(MappedRead const * const read, int reason) {
	AtomicInc(&m_UnmappedReads);
	Log.Debug(LOG_OUTPUT_DETAILS, "Read %s (%i) not mapped (%i)", read->name, read->ReadId, reason);
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
		m_CurStart = 0;
		NGMSignal(&m_CSWait);
	}

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

NGMTHREADFUNC _NGM::ThreadFunc(void* data) {
	NGMTask * task = (NGMTask*) data;
	int tid = task->m_TID;

	try {
		Log.Verbose("Running thread %i", tid);

		task->Run();

		Log.Verbose("Thread %i run return, finishing", tid);

		NGM.FinishThread(tid);

		Log.Verbose("Thread %i finished", tid);
	} catch (...) {
		Log.Error("Unhandled exception in thread %i", tid);
	}

	Log.Verbose("ThreadFunc on thread %i returning", tid);

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

void _NGM::StartThreads() {
	int threadcount = Config.GetInt("cpu_threads", 1, 0);

#ifdef _DEBUGCS
	threadcount = 1;
#endif

	StartCS(threadcount);

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

IAlignment * _NGM::CreateAlignment(int const mode) {
	IAlignment * instance = 0;

	if (Config.GetInt("affine")) {
		if (Config.Exists("gpu")) {
			Log.Error("GPU doesn't support affine gap penalties. Please remove --affine or -g/--gpu.");
			Fatal();
		}
		instance = new EndToEndAffine();
	} else {
		int dev_type = CL_DEVICE_TYPE_CPU;

		if (Config.Exists("gpu")) {
			dev_type = CL_DEVICE_TYPE_GPU;
		}

		Log.Debug(LOG_INFO, "INFO", "Mode: %d GPU: %d", mode, mode & 0xFF);
		OclHost * host = new OclHost(dev_type, mode & 0xFF, Config.GetInt("cpu_threads"));

		int ReportType = (mode >> 8) & 0xFF;
		switch (ReportType) {
			case 1:

			instance = new SWOclCigar(host);
			break;
			default:
			Log.Error("Unsupported report type %i", mode);
			break;
		}
	}

	return instance;
}

void _NGM::DeleteAlignment(IAlignment* instance) {
	OclHost * host = 0;
	if (!Config.GetInt("affine")) {
		host = ((SWOcl *) instance)->getHost();
	}
	Log.Verbose("Delete alignment called");
	if (instance != 0) {
		delete instance;
		instance = 0;
	}

	if (host != 0) {
		delete host;
		host = 0;
	}
}

void _NGM::MainLoop() {

	bool const isPaired = Config.GetInt("paired") > 0;
	int const threadcount = Config.GetInt("cpu_threads", 1, 0);
	bool const progress = Config.GetInt("no_progress") != 1;
	bool const argos = Config.Exists(ARGOS);
	long lastProcessd = 0;
	while (Running()) {
		Sleep(50);
		if (progress) {
			int processed = std::max(1, NGM.GetMappedReadCount() + NGM.GetUnmappedReadCount());
			if(!argos || (lastProcessd + 1000000) < processed) {
				lastProcessd = processed;
				if (!isPaired) {
					Log.Progress("Mapped: %d, CMR/R: %d, CS: %d (%d), R/S: %d, Time: %.2f %.2f %.2f", processed, NGM.Stats->avgnCRMS, NGM.Stats->csLength, NGM.Stats->csOverflows, NGM.Stats->readsPerSecond * threadcount, NGM.Stats->csTime, NGM.Stats->scoreTime, NGM.Stats->alignTime);
				} else {
					Log.Progress("Mapped: %d, CMR/R: %d, CS: %d (%d), R/S: %d, Time: %.2f %.2f %.2f, Pairs: %.2f %.2f", processed, NGM.Stats->avgnCRMS, NGM.Stats->csLength, NGM.Stats->csOverflows, NGM.Stats->readsPerSecond * threadcount, NGM.Stats->csTime, NGM.Stats->scoreTime, NGM.Stats->alignTime, NGM.Stats->validPairs, NGM.Stats->insertSize);
				}
			}
		}
	}
}

