#include "Aligner.h"
#include "NGM.h"

#undef module_name
#define module_name "ALIGN"

//#ifdef STATIC

#include "OclHost.h"
#include "SWOclCigar.h"

IAlignment * CreateAlignment(int const mode) {
	//int dev_type = Config.GetInt("ocl_device");
	int dev_type = CL_DEVICE_TYPE_GPU;
#ifdef CPU
	dev_type = CL_DEVICE_TYPE_CPU;
#else
#endif
#ifndef NDEBUG
	Log.Error("Mode: %d GPU: %d", mode, mode & 0xFF);
#endif
	OclHost * host = new OclHost(dev_type, mode & 0xFF, Config.GetInt("ocl_threads"));

	SWOcl * instance = 0;

//#ifndef NDEBUG
	//Log.Error("Alignment mode: %d", mode);
//#endif
	int ReportType = (mode >> 8) & 0xFF;
	switch (ReportType) {
	case 0:
#ifndef NDEBUG
		Log.Message("Output: text");
#endif
//		instance = new SWOclAlignment(host);
		//			instance = new SWOclCigar(host);
		break;
	case 1:
#ifndef NDEBUG
		Log.Message("Output: cigar");
#endif
		instance = new SWOclCigar(host);
		break;
	default:
		Log.Error("Unsupported report type %i", mode);
		break;
	}
	return instance;
}

void DeleteAlignment(IAlignment* instance) {
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
		delete host;
		host = 0;
	}
}
//#endif

/*
 *  Aligner ist eine Thread/CUDA/GPU-Einheit zur Berechnung von Scores/Alignments
 *
 *
 */
class AlignmentDispatcher::Aligner {
private:
	NGMMutex m_Mutex;
	NGMThreadWait m_EntryWait;
	NGMThreadWait m_ReturnWait;

	int const m_GPU;
	int const m_OutputFormat;

	IAlignment * m_Kernel;
	volatile bool m_KernelInit;

	bool m_Running;
	bool m_Finished;

	int m_IvkType;

	int m_IvkMode;
	int m_IvkBatchSize;
	char const * const * m_IvkRefSeqList;
	char const * const * m_IvkQrySeqList;
	void * m_IvkResults;
	void * m_IvkData;
	int m_IvkReturn;

	static NGMTHREADFUNC Run(void* pObj) {
		Aligner * align = (Aligner*) pObj;

		Log.Message("Alignment Thread for GPU %i startet", align->m_GPU);
		NGMLock(&align->m_Mutex);
		Log.Verbose("Creating alignment kernel on GPU %i", align->m_GPU);
		align->m_Kernel = CreateAlignment(align->m_GPU | (std::min(align->m_OutputFormat, 1) << 8));
		Log.Verbose("Alignment kernel on GPU %i created", align->m_GPU);
		if (align->m_Kernel == 0)
			Fatal();

		Log.Message("Alignment created for GPU %i.", align->m_GPU);
		NGMSignal(&align->m_ReturnWait);

		while (true) {
			while (align->m_IvkType <= 0) {
				if (align->m_Finished) {
					NGMUnlock(&align->m_Mutex);
					Log.Message("Alignment Thread for GPU %i finished", align->m_GPU);
					return 0;
				}

				Log.Verbose("Alignment on GPU %i waiting for jobs...", align->m_GPU);
				NGMWait(&align->m_Mutex, &align->m_EntryWait);
			}
			Log.Verbose("Launching Kernel on GPU %i...", align->m_GPU);

//			for(int i = 0; i < align->m_IvkBatchSize; ++i) {
//				Log.Message("%c", ((char *)align->m_IvkData)[i] + '0');
//			}

			if (align->m_IvkType == 1)
				align->m_IvkReturn = align->m_Kernel->BatchScore(align->m_IvkMode, align->m_IvkBatchSize, align->m_IvkRefSeqList,
						align->m_IvkQrySeqList, (float*) align->m_IvkResults, align->m_IvkData);
			else if (align->m_IvkType == 2)
				align->m_IvkReturn = align->m_Kernel->BatchAlign(align->m_IvkMode, align->m_IvkBatchSize, align->m_IvkRefSeqList,
						align->m_IvkQrySeqList, (Align*) align->m_IvkResults, align->m_IvkData);
			else
				align->m_IvkReturn = 0;

			align->m_IvkType = 0;

			Log.Verbose("Kernel returned (%i) on GPU %i", align->m_IvkReturn, align->m_GPU);
			NGMSignal(&align->m_ReturnWait);
			Log.Verbose("signalled");

		}
	}
public:
	bool InUse;
	ulong BusyTime;

	//static pfCreateAlignment CreateAlignment;
	//static pfDeleteAlignment DeleteAlignment;

	Aligner(int const gpu) :
			m_Mutex(), m_EntryWait(), m_ReturnWait(), m_GPU(gpu), m_OutputFormat(Config.GetInt("format", 0, 2)), m_Kernel(0), m_KernelInit(
					false), m_Running(false), m_Finished(false), m_IvkType(0), InUse(false), BusyTime(0) {
		NGMInitMutex(&m_Mutex);
		NGMInitWait(&m_EntryWait);
		NGMInitWait(&m_ReturnWait);

		NGMCreateThread(&Aligner::Run, this, true);
	}

	int Dispatch(int type, int const mode, int const batchSize, char const * const * const refSeqList,
			char const * const * const qrySeqList, void * const results, void * extData) {
		NGMLock(&m_Mutex);
		if (m_Running) {
			NGMUnlock(&m_Mutex);
			return -1;
		}

		m_IvkType = type;
		m_IvkMode = mode | (std::min(m_OutputFormat, 1) << 8);
		m_IvkBatchSize = batchSize;
		m_IvkRefSeqList = refSeqList;
		m_IvkQrySeqList = qrySeqList;
		m_IvkResults = results;
		m_IvkData = extData;
		m_IvkReturn = 0;

		m_Running = true;

		NGMSignal(&m_EntryWait);

		while (m_IvkType > 0) {
			//Log.Error("Waiting");
			NGMWait(&m_Mutex, &m_ReturnWait);
		}

		int ret = m_IvkReturn;
		m_Running = false;
		NGMUnlock(&m_Mutex);

		return ret;
	}

	IAlignment const * Kernel() {
		if (m_Kernel == 0) {
			NGMLock(&m_Mutex);
			while (m_Kernel == 0)
				NGMWait(&m_Mutex, &m_ReturnWait);
			NGMUnlock(&m_Mutex);
		}
		return m_Kernel;
	}

	void Shutdown() {
		NGMLock(&m_Mutex);
		m_Finished = true;
		NGMSignal(&m_EntryWait);
		NGMUnlock(&m_Mutex);
	}
};

//pfCreateAlignment AlignmentDispatcher::Aligner::CreateAlignment = 0;
//pfDeleteAlignment AlignmentDispatcher::Aligner::DeleteAlignment = 0;

AlignmentDispatcher * AlignmentDispatcher::sInstance = 0;

NGMOnceControl dispatcher_once = NGM_ONCE_INIT;

void AlignmentDispatcher::CreateInstance() {
	sInstance = new AlignmentDispatcher();
}

AlignmentDispatcher * AlignmentDispatcher::Instance() {
	if (sInstance == 0) {
		NGMOnce(&dispatcher_once, &AlignmentDispatcher::CreateInstance);
	}
	return sInstance;
}
void AlignmentDispatcher::Shutdown() {
	for (int i = 0; i < m_TotalAligner; ++i) {
		m_Aligner[i]->Shutdown();
		delete m_Aligner[i];
	}
	delete this;
	sInstance = 0;
}

AlignmentDispatcher::AlignmentDispatcher() :
		m_SelectionMutex(), m_SelectionWait(), m_CPUEnabled(false), m_CPUKernels(), m_CPUSelectionMutex(), m_StartTime(GetTickCount()) {
	NGMInitMutex(&m_SelectionMutex);
	NGMInitWait(&m_SelectionWait);
	NGMInitMutex(&m_CPUSelectionMutex);

	m_Aligner[0] = 0;

	m_TotalAligner = m_FreeAligner = InitAligners();

	if (m_TotalAligner == 0 && !m_CPUEnabled) {
		Log.Error("Both CPU and GPU DLL missing or failed to load.");
		Fatal();
	} else {
		if (m_TotalAligner > 0)
			Log.Message("Dispatcher with %i GPU Aligners%s created", m_TotalAligner, (m_CPUEnabled) ? " and additional CPU Aligners" : "");
		else
			Log.Message("Dispatcher with CPU Aligners created");
	}
}

AlignmentDispatcher::~AlignmentDispatcher() {
}

int AlignmentDispatcher::InitAligners() {
	typedef int (*pfCookie)();
	int threadcount = 0;

	/*** Load GPU DLL ***/
	if (Config.Exists("gpu") && Config.GetInt("gpu", 0, cMaxAligner) > 0) {
		char const * gpuDlls[2] = { "libMASonOpenClGpu.so", "libMASonCUDA.so" };
		char const * gpuDll = 0;
		if (!Config.Exists("gpu_dll")) {
			if (Config.GetInt("bs_mapping") == 0) {
				Log.Message("Detecting alignment dll.");
				if (CheckDLL(gpuDlls[0])) {
					gpuDll = gpuDlls[0];
					Log.Message("Using %s", gpuDll);
				} else if (CheckDLL(gpuDlls[1])) {
					gpuDll = gpuDlls[1];
					Log.Message("Using %s", gpuDll);
				} else {
					Log.Error("Error loading GPU DLL. If the MASon libs or not located in standard directories, pleas specifiy mason_path in the config file.");
					Fatal();
				}
			} else {
				Log.Message("Detecting alignment dll for bs-mapping.");
				if (CheckDLL(gpuDlls[1])) {
					gpuDll = gpuDlls[1];
					Log.Message("Using %s", gpuDll);
				} else {
					Log.Error("Error loading GPU DLL. If the MASon libs or not located in standard directories, pleas specifiy mason_path in the config file.");
					Fatal();
				}
			}
		} else {
			gpuDll = Config.GetString("gpu_dll");
			Log.Message("GPU dll specified. Using %s", gpuDll);
		}
		int const dll = InitDLL(gpuDll);
		if (dll == -1) {
			Log.Error("Error loading GPU DLL");
			Fatal();
		} else {
			pfCookie pfC = (pfCookie) GetDLLFunc(dll, "Cookie", true);
			if (pfC() != cCookie) {
				Log.Error("Invalid DLL version. Please update and rebuild the DLL from SVN");
			} else {
				//Aligner::CreateAlignment = (pfCreateAlignment) GetDLLFunc(dll, "CreateAlignment", true);
				//Aligner::DeleteAlignment = (pfDeleteAlignment) GetDLLFunc(dll, "DeleteAlignment", false);

				threadcount = Config.GetInt("gpu", 1, cMaxAligner);
				Log.Message("Alignment Threads: %i", threadcount);
				int * gpus = new int[threadcount];
				Config.GetIntArray("gpu", gpus, threadcount);

				for (int i = 0; i < threadcount; ++i) {
					Log.Green("Using GPU %d for alignment/score computation.", gpus[i]);
					m_Aligner[i] = new Aligner(gpus[i]);
				}

				delete[] gpus;

				m_ScoreBatchSize = m_Aligner[0]->Kernel()->GetScoreBatchSize();
				m_AlignBatchSize = m_Aligner[0]->Kernel()->GetAlignBatchSize();
			}
		}

		/*** Load CPU DLL ***/
		//TODO: detect dlls
	} else {
		char const * cpuDlls[2] = { "libMASonOpenClCpu.so", "libMASonSSE.so" };
		char const * cpuDll = 0;
		if (!Config.Exists("cpu_dll")) {
			Log.Message("Detecting alignment dll.");
			if (CheckDLL(cpuDlls[0])) {
				cpuDll = cpuDlls[0];
				Log.Message("Using %s", cpuDll);
			} else if (CheckDLL(cpuDlls[1])) {
				cpuDll = cpuDlls[1];
				Log.Message("Using %s", cpuDll);
			} else {
				Log.Error("Error loading CPU DLL. If the MASon libs or not located in standard directories, pleas specifiy mason_path in the config file.");
				Fatal();
			}
		} else {
			cpuDll = Config.GetString("cpu_dll");
			Log.Message("CPU dll specified. Using %s", cpuDll);
		}
		int const cpu_dll = InitDLL(cpuDll);
		if (cpu_dll == -1)
			Log.Error("Error loading CPU DLL");
		else {
			pfCookie cpuCookie = (pfCookie) GetDLLFunc(cpu_dll, "Cookie", true);
			if (cpuCookie() != cCookie) {
				Log.Error("Invalid CPU DLL version. Please update and rebuild the DLL from SVN");
			} else {
				m_CreateCPUAlignment = (pfCreateAlignment) GetDLLFunc(cpu_dll, "CreateAlignment", true);
				m_DeleteCPUAlignment = (pfDeleteAlignment) GetDLLFunc(cpu_dll, "DeleteAlignment", false);
				m_CPUEnabled = true;

				if (threadcount == 0) {
					IAlignment * aligner = GetCPUAligner();
					m_ScoreBatchSize = aligner->GetScoreBatchSize();
					m_AlignBatchSize = aligner->GetAlignBatchSize();
				}
			}
		}
	}
	/*********************/

	return threadcount;
}

int AlignmentDispatcher::GetScoreBatchSize() const {
	return m_ScoreBatchSize;
}

int AlignmentDispatcher::GetAlignBatchSize() const {
	return m_AlignBatchSize;
}

int AlignmentDispatcher::BatchScore(int const mode, int const batchSize, char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results, void * extData) {
	Aligner * aligner = PickFreeAligner(0);
	if (aligner == 0 && extData != 0)
		Log.Error("Warning: CPU Calculation with additional data (extData) is not supported and may yield wrong results");

	int bSize = batchSize;
	static int sBatchLimit = GetScoreBatchSize();
	if (bSize > sBatchLimit)
		bSize = sBatchLimit;

	int res = 0;
	int cur = 0;

	while (aligner == 0 && bSize > 0) {
		IAlignment * cpu = GetCPUAligner();
		int cpu_batchsize = cpu->GetScoreBatchSize();
		if (cpu_batchsize > bSize)
			cpu_batchsize = bSize;

		res += cpu->BatchScore(mode, cpu_batchsize, refSeqList + cur, qrySeqList + cur, results + cur, extData);
		cur += cpu_batchsize;
		bSize -= cpu_batchsize;

		FreeCPUAligner(cpu);
		aligner = PickFreeAligner(cur);
	}

	if (bSize > 0) {
		ulong startTime = GetTickCount();
		res += aligner->Dispatch(1, mode, bSize, refSeqList, qrySeqList, results, extData);
		ulong deltaT = startTime - GetTickCount();
		aligner->BusyTime += deltaT;

		aligner->InUse = false;
		AtomicInc(&m_FreeAligner);
		NGMSignal(&m_SelectionWait);
	}
	return res;
}

int AlignmentDispatcher::BatchAlign(int const mode, int const batchSize, char const * const * const refSeqList,
		char const * const * const qrySeqList, Align * const results, void * extData) {
	Aligner * aligner = PickFreeAligner(0);
	if (aligner == 0 && extData != 0)
		Log.Error("Warning: CPU Calculation with additional data (extData) is not supported and may yield wrong results");

	static int sBatchLimit = GetAlignBatchSize();
	int bSize = batchSize;
	if (bSize > sBatchLimit)
		bSize = sBatchLimit;

	int res = 0;
	int cur = 0;

	while (aligner == 0 && bSize > 0) {
		IAlignment * cpu = GetCPUAligner();
		int cpu_batchsize = cpu->GetScoreBatchSize();
		if (cpu_batchsize > bSize)
			cpu_batchsize = bSize;

		res += cpu->BatchAlign(mode, cpu_batchsize, refSeqList + cur, qrySeqList + cur, results + cur, extData);
		cur += cpu_batchsize;
		bSize -= cpu_batchsize;

		FreeCPUAligner(cpu);
		if (bSize > 0)
			aligner = PickFreeAligner(cur);
	}

	if (bSize > 0) {
		ulong startTime = GetTickCount();
		res += aligner->Dispatch(2, mode, bSize, refSeqList, qrySeqList, results, extData);
		ulong deltaT = startTime - GetTickCount();
		aligner->BusyTime += deltaT;

		aligner->InUse = false;
		AtomicInc(&m_FreeAligner);
		NGMSignal(&m_SelectionWait);
	}
	return res;
}

bool AlignmentDispatcher::ChooseCPU(int desire) {
	// Keine GPU-Kernel vorhanden - immer CPU wählen
	if (m_TotalAligner == 0)
		return true;

	// Keine CPU-Kernel vorhanden - immer GPU wählen
	if (!m_CPUEnabled)
		return false;

	// GPU frei - diese wählen
	if (m_FreeAligner > 0)
		return false;

	if (desire == 0) // Voller GPU-Batch
			{
		if (GPULoad() > 0.8f)
			return true;
	} else if (desire > 0) {
		return true;
	}
	return false;
}

AlignmentDispatcher::Aligner * AlignmentDispatcher::PickFreeAligner(int desire) {
	Aligner * aligner = 0;

	if (ChooseCPU(desire))
		return 0;

	// Pick free GPU Kernel
	NGMLock(&m_SelectionMutex);
	while (!(m_FreeAligner > 0))
		NGMWait(&m_SelectionMutex, &m_SelectionWait);

	int i = 0;
	while (aligner == 0 && i < m_TotalAligner) {
		if (!m_Aligner[i]->InUse)
			aligner = m_Aligner[i];
		else
			++i;
	}
	if (aligner == 0)
		Log.Error("No aligner found. Crashing...");

	aligner->InUse = true;
	AtomicDec(&m_FreeAligner);
	NGMUnlock(&m_SelectionMutex);
	return aligner;
}

// Resets GPU load calculculations
void AlignmentDispatcher::ResetLoad() {
	for (int i = 0; i < m_TotalAligner; ++i) {
		m_Aligner[i]->BusyTime = 0;
	}
	m_StartTime = GetTickCount();
}

// Returns the average load of all GPUs in Range 0f..1f (0-100%)
float AlignmentDispatcher::GPULoad() const {
	ulong totalGPUTime = 0;
	for (int i = 0; i < m_TotalAligner; ++i) {
		totalGPUTime += m_Aligner[i]->BusyTime;
	}

	float load = ((float) totalGPUTime / (float) m_TotalAligner) / (float) (GetTickCount() - m_StartTime);
	if (load > 1.0f)
		load = 1.0f;
	if (load < 0)
		load = 0;

	return load;
}

// Picks a free CPU Kernel or, if none are available, creates one
IAlignment * AlignmentDispatcher::GetCPUAligner() {
	IAlignment * aligner = 0;
	static int i = 0;
	if (!m_CPUEnabled)
		Log.Error("Tried to get CPU Kernel without CPU DLL.");

	NGMLock(&m_CPUSelectionMutex);
	if (m_CPUKernels.size() > 0) {
		aligner = m_CPUKernels[m_CPUKernels.size() - 1];
		m_CPUKernels.pop_back();
	} else {
		aligner = m_CreateCPUAlignment(i++ | (std::min(Config.GetInt("format", 0, 2), 1) << 8));
	}
	NGMUnlock(&m_CPUSelectionMutex);
	return aligner;
}

void AlignmentDispatcher::FreeCPUAligner(IAlignment * aligner) {
	NGMLock(&m_CPUSelectionMutex);
	m_CPUKernels.push_back(aligner);
	NGMUnlock(&m_CPUSelectionMutex);
}
