/*
 * SWOcl_export.cpp
 *
 *  Created on: Jun 18, 2011
 *      Author: philipp_
 */

//#define _TEST
//#define _DEBUG

#include "OclHost.h"
#include "SWOclCigar.h"
#include "SWOclAlignment.h"

#ifdef _WIN32
#define dllexport  __declspec(dllexport)
#else
#define dllexport
#endif

// Logging:
ILog const * _log = 0;
extern "C" dllexport void SetLog(ILog * log) {
	_log = log;
}

extern "C" dllexport int Cookie() {
	return 0x10201130;
}

IConfig * _config = 0;
extern "C" dllexport void SetConfig(IConfig * config) {
	_config = config;
}

extern "C" dllexport bool IsAvailable() {
	return true;
}

extern "C" dllexport IAlignment * CreateAlignment(int const mode) {
	//int dev_type = Config.GetInt("ocl_device");
	int dev_type = CL_DEVICE_TYPE_GPU;
#ifdef CPU
	dev_type = CL_DEVICE_TYPE_CPU;
#else
#endif
#ifndef NDEBUG
	Log.Error("Mode: %d GPU: %d", mode, mode & 0xFF);
#endif
	OclHost * host = new OclHost(dev_type, mode & 0xFF, Config.GetInt(
			"ocl_threads"));

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
		instance = new SWOclAlignment(host);
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
	}
	return instance;
}

extern "C" dllexport void ExternalDeleteString(char* mem) {
	delete[] mem;
}

extern "C" dllexport void DeleteAlignment(SWOcl* instance) {
	OclHost * host = instance->getHost();
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

