/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "NGM.h"

#include <stdlib.h>
#include <stdio.h>
#ifdef _WIN32
#include <conio.h>
#endif
#ifdef _WIN32
#include <crtdbg.h>
#endif
#include <unistd.h>
#include <sstream>

#include "Log.h"
#include "CS.h"
#include "Version.h"
#include "Timing.h"
//#include "UpdateCheck.h"
#include "ArgParser.h"

#undef module_name
#define module_name "MAIN"

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

#ifdef NDEBUG
bool cDebug = false;
#else
bool cDebug = true;
#endif

ILog const * _log = 0;
IConfig * _config = 0;


int main(int argc, char * argv[]) {

	std::stringstream version;
	version << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_BUILD;

	Log.Message("ngmlr %s%s (build: %s %s, start: %s)", version.str().c_str(), cDebug ? " debug" : "", __DATE__, __TIME__, currentDateTime().c_str());
	Log.Message("Contact: fritz.sedlazeck@gmail.com, philipp.rescheneder@gmail.com");

	Timer tmr;
	tmr.ST();

	_NGM::AppName = argv[0];

	_config = new ArgParser(argc, argv); // Parses command line & parameter file
	_log = &Log;

#ifdef DEBUGLOG
	if (Config.getLog()) {
		Log.Message("Writing debug log to stdout. Please use -o/--output for SAM/BAM output.");
		//Init checks if first parameter is != 0. Thus "LOG" is passed as a dummy string.
		_Log::Init("LOG", Config.getLogLevel());// Inits logging to file
	}
#else
	_Log::Init(0, 0); // Inits logging to file
#endif

//	if (Config.getUpdateCheck()) {
//		UpdateCheckInterface::remoteCheck();
//		exit(0);
//	}

	_Log::setColor(Config.getColor());


	NGM; // Init Core

	NGM.InitProviders();

	try {
		NGM.StartThreads();

		bool const progress = Config.getProgress();

		if(!progress) {
			Log.Message("No progress information (use --progress)");
		}
		Log.Message("Mapping reads...");
		if(Config.isStdIn()) {
			Log.Message("Waiting for data from stdin");
		}


		NGM.MainLoop();

		int discarded = NGM.GetReadReadCount() - (NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount());
		if (discarded != 0) {
			Log.Message("Done (%i reads mapped (%.2f%%), %i reads not mapped (%i discarded), %i lines written)(elapsed: %dm, %d r/s)", NGM.GetMappedReadCount(), (float)NGM.GetMappedReadCount() * 100.0f / (float)(std::max(1, NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount() + discarded)),NGM.GetUnmappedReadCount() + discarded, discarded, NGM.GetWrittenReadCount(), (int) (tmr.ET() / 60.0f), (int)(NGM.GetMappedReadCount() * 1.0f / tmr.ET()));
		} else {
			Log.Message("Done (%i reads mapped (%.2f%%), %i reads not mapped, %i lines written)(elapsed: %dm, %d r/s)", NGM.GetMappedReadCount(), (float)NGM.GetMappedReadCount() * 100.0f / (float)(std::max(1, NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount())),NGM.GetUnmappedReadCount(),NGM.GetWrittenReadCount(), (int) (tmr.ET() / 60.0f), (int)(NGM.GetMappedReadCount() * 1.0f / tmr.ET()));
		}

		NGM.StopThreads();

	} catch (...) {
		Log.Error("Unhandled exception in control thread.");
	}



//	if (! Config.getUpdateCheck()) {
//		UpdateCheckInterface::reminder();
//	}

	CS::Cleanup();
	_SequenceProvider::Cleanup();
	delete _config;
	delete _NGM::pInstance;
	_Log::Cleanup();

	CleanupPlatform();

	return 0;
}


// actually platform specific.../care
uloc const FileSize(char const * const filename) {
	FILE * fp = fopen(filename, "rb");
	if (fp == 0) {
		Log.Warning("Tried to get size of nonexistent file %s", filename);
		return 0;
	}

	if (fseek(fp, 0, SEEK_END) != 0)
		return 0;

#ifdef __APPLE__
	uloc end = ftello(fp);
#else
	uloc end = ftello64(fp);
#endif

	fclose(fp);
	return end;
}
