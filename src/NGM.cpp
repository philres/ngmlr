/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "NGM.h"

#include <memory.h>
#include <stdlib.h>
#include <limits.h>

#include "CS.h"
#include "PrefixTable.h"
#include "ReadProvider.h"
#include "PlainFileWriter.h"
#include "SAMWriter.h"
#include "Timing.h"
#include "StrippedSW.h"

#undef module_name
#define module_name "NGM"

_NGM * _NGM::pInstance = 0;
NGMOnceControl _NGM::once_control = NGM_ONCE_INIT;
char const * _NGM::AppName = 0;

namespace __NGM {
	inline int min(int a, int b) {
		return (a < b) ? a : b;
	}
}

void _NGM::Init() {
	pInstance = new _NGM();
}

_NGM & _NGM::Instance() {
	NGMOnce(&_NGM::once_control, Init);
	return *pInstance;
}

_NGM::_NGM() : Stats(new NGMStats()), m_ActiveThreads(0), m_NextThread(0), m_CurStart(0), m_CurCount(0), m_SchedulerMutex(), m_SchedulerWait(), m_TrackUnmappedReads(false), m_UnmappedReads(0), m_MappedReads(0), m_WrittenReads(0), m_ReadReads(0), m_RefProvider(0), m_ReadProvider(0) {

	char const * const output_name = Config.getOutputFile();
//	if (!Config.getBAM()) {
	if (output_name != 0) {
		Log.Message("Opening for output (SAM): %s", output_name);
	} else {
		Log.Message("Writing output (SAM) to stdout");
	}

	m_Output = new PlainFileWriter(output_name);

//	} else {
//		if (output_name != 0) {
//			Log.Message("Opening for output (BAM): %s", output_name);
//		} else {
//			Log.Message("Wrinting output (BAM) to stdout");
//		}
//		m_Output = new FileWriterBam(output_name);
//	}

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

}

void _NGM::InitProviders() {
	CS::Init();

	SequenceProvider.Init(); // Prepares input data

	m_RefProvider = new CompactPrefixTable();

	m_ReadProvider = new ReadProvider();
	uint readCount = m_ReadProvider->init();
}

_NGM::~_NGM() {

	delete Stats;
	Stats = 0;

	if (m_RefProvider != 0)
		delete m_RefProvider;

	if (m_ReadProvider != 0)
		delete m_ReadProvider;
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
//		if (!Config.getBAM()) {
		delete (FileWriter*) m_Output;
//		} else {
//			delete (FileWriterBam*) m_Output;
//		}
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


	desBatchSize &= ~1;

	if (m_CurCount == 0) {
		m_CurStart = 0;
		NGMSignal(&m_CSWait);
	}

	list.reserve(desBatchSize);
	int count = 0;

	//Long PacBio reads are split into smaller parts.
	//Each part should have its own id.
	int idJump = 2000;

	int i = 0;
	while (count < desBatchSize && !eof) {
		MappedRead * read1 = 0;
		eof = !NGM.GetReadProvider()->GenerateRead(m_CurStart + i * idJump, read1, 0, read1);
		i += 1;
		if (!eof) {
			if (read1 != 0) {
				count += 1;
				Stats->readLengthSum += read1->length;
				if (read1->group == 0) {
					// Short read found: not split into read group
					list.push_back(read1);
				} else {
					// Long read found: push subreads

					for (int j = 0; j < read1->group->readNumber; ++j) {
						list.push_back(read1->group->reads[j]);
					}
				}
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
		NGM.FinishThread(tid);
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

void _NGM::StopThreads() {

}

void _NGM::StartThreads() {

	StartCS(Config.getThreads());

}

void _NGM::StartCS(int cs_threadcount) {

	for (int i = 0; i < cs_threadcount; ++i) {
		NGMTask * cs = 0;
		cs = new CS();
		StartThread(cs, -1);
	}

}

IAlignment * _NGM::CreateAlignment(int const mode) {
	IAlignment * instance = 0;

	switch (Config.getSubreadAligner()) {
	case 2:
		instance = new StrippedSW();
		break;
	default:
		Log.Error("Invalid subread alignerd: %d", Config.getSubreadAligner());
		throw "";
	}
	return instance;
}

void _NGM::DeleteAlignment(IAlignment* instance) {

	Log.Verbose("Delete alignment called");
	if (instance != 0) {
		delete instance;
		instance = 0;
	}

}

void _NGM::MainLoop() {

	Timer tmr;
	tmr.ST();

	bool const progress = Config.getProgress();

	int processed = 0;
	float runTime = 0.0f;
	float readsPerSecond = 0.0f;
	float alignSuccessRatio = 0.0f;
	int avgCorridor = 0;

	while (Running()) {
		Sleep(2000);
		if (progress) {
			processed = std::max(1, NGM.GetMappedReadCount() + NGM.GetUnmappedReadCount());

			runTime = tmr.ET();
			readsPerSecond = processed / runTime;
			if((NGM.Stats->alignmentCount + NGM.Stats->invalidAligmentCount) > 0) {
				avgCorridor = NGM.Stats->corridorLen / (NGM.Stats->alignmentCount + NGM.Stats->invalidAligmentCount);
				alignSuccessRatio = NGM.Stats->alignmentCount * 1.0f / (NGM.Stats->alignmentCount + NGM.Stats->invalidAligmentCount);
			}
			float avgAlignPerc = 0.0f;
			int avgReadLenght = 0;
			float alignRate = 0.0f;
			if(processed  > 0) {
				avgAlignPerc = NGM.Stats->avgAlignPerc / std::max(1, NGM.GetMappedReadCount());
				avgReadLenght = (int)(NGM.Stats->readLengthSum / processed);
				alignRate = NGM.GetMappedReadCount() * 1.0f / processed;
			}
			Log.Progress("Processed: %d (%.2f), R/S: %.2f, RL: %d, Time: %.2f %.2f %.2f, Align: %.2f, %d, %.2f", processed, alignRate, readsPerSecond, avgReadLenght, NGM.Stats->csTime, NGM.Stats->scoreTime, NGM.Stats->alignTime, alignSuccessRatio, avgCorridor, avgAlignPerc);
		}
	}

	if (progress) {
		processed = std::max(1, NGM.GetMappedReadCount() + NGM.GetUnmappedReadCount());
		runTime = tmr.ET();
		readsPerSecond = processed / runTime;
		if((NGM.Stats->alignmentCount + NGM.Stats->invalidAligmentCount) > 0) {
			alignSuccessRatio = NGM.Stats->alignmentCount * 1.0f / (NGM.Stats->alignmentCount + NGM.Stats->invalidAligmentCount);
			avgCorridor = NGM.Stats->corridorLen / (NGM.Stats->alignmentCount + NGM.Stats->invalidAligmentCount);
		}
		float avgAlignPerc = 0.0f;
		int avgReadLenght = 0;
		float alignRate = 0.0f;
		if(processed > 0) {
			avgAlignPerc = NGM.Stats->avgAlignPerc / std::max(1, NGM.GetMappedReadCount());
			avgReadLenght = (int)(NGM.Stats->readLengthSum / processed);
			alignRate = NGM.GetMappedReadCount() * 1.0f / processed;
		}
		Log.Message("Processed: %d (%.2f), R/S: %.2f, RL: %d, Time: %.2f %.2f %.2f, Align: %.2f, %d, %.2f", processed, alignRate, readsPerSecond, avgReadLenght, NGM.Stats->csTime, NGM.Stats->scoreTime, NGM.Stats->alignTime, alignSuccessRatio, avgCorridor, avgAlignPerc);
	}
}

