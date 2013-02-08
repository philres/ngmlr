//#include "NGM.h"

#include "Log.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "PlatformSpecifics.h"
//#include "Config.h"

#include "NGMThreads.h"

#include <vector>
#include <string>

#include <ctime>

#include "Debug.h"

namespace __Log {
_Log * pInstance = 0;
NGMMutex mutex;
NGMOnceControl once_control = NGM_ONCE_INIT;

int rwd;

int warningCount = 0;

void Init() {
	pInstance = new _Log();
	InitConsole();

	NGMInitMutex(&mutex);
}

int filterlvl = 0;
bool init = false;
char preBuffer[1024];
bool preInit = true;
bool logToFile = false;
std::vector<std::string> & msgLog() {
	static std::vector<std::string> * res = new std::vector<std::string>();
	return *res;
}

FILE * fp = 0;

char const * lvlStr[] = { "", "[WARNING]", "[ERROR]", "" };

void LogToFile(int lvl, char const * const title, char const * const s, va_list args) {
	if (!logToFile && !preInit)
		return;

	int written = 0;

	if (title != 0)
		written += sprintf(preBuffer, "%s[%s] ", lvlStr[lvl], title);
	else
		written += sprintf(preBuffer, "%s ", lvlStr[lvl]);

	vsprintf(preBuffer + written, s, args);

	if (preInit)
		msgLog().push_back(std::string(preBuffer));

	if (logToFile) {
		fprintf(fp, "%s\n", preBuffer);
		fflush(fp);
	}
}

// rewind lines
inline void rwl() {
	for (int i = 0; i < rwd; ++i)
		printf("\033[A\033[2K");
	rwd = 0;
}

void LogToConsole(int lvl, char const * const title, char const * const s, va_list args) {
	if (lvl < filterlvl)
		return;

	bool progress = lvl == 99;
	if (progress)
		lvl = 0;
	rwl();

	if (title != 0) {
		if (lvl > 0)
			SetConsoleColor((ConsoleColor) (MessageTitle + (lvl * 2)));
		printf("[%s] ", title);
	}

	SetConsoleColor((ConsoleColor) (Message + (lvl * 2)));
	if (args != 0)
		vprintf(s, args);
	else
		printf("%s", s);
	ResetConsoleColor();
	printf("\n");
	if (progress) {
		rwd = 1;
	}
}

}

using namespace __Log;

std::string add_timestamp(std::string str) {
	std::string::size_type pos = str.find("%s");
	if (pos != std::string::npos) {
		//char months[] = "Jan\0Feb\0Mar\0Apr\0May\0Jun\0Jul\0Aug\0Sep\0Oct\0Nov\0Dec\0";
		const time_t t = time(0);
		tm * tm = localtime(&t);

		char timestamp[20];
		sprintf(timestamp, "%4i-%02i-%02i_%02i-%02i-%02i", tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday, tm->tm_hour, tm->tm_min,
				tm->tm_sec); //months+m*4

		return str.replace(pos, 2, timestamp);
	}
	return str;
}

void _Log::Init() {
	return; // File logging disbaled

	init = true;
	NGMLock(&__Log::mutex);

	try {
//		if (Config.Exists("log")) {
//			if ((logToFile = (Config.GetInt("log") != 0))) {
//				char const * filename = 0;
//				if (Config.Exists("logfile"))
//					filename = Config.GetString("logfile");
//
//				if (filename == 0) {
//					filename = "ngm_%s.log";
//					filename = add_timestamp(filename).c_str();
//				}
//
//				fp = fopen(filename, "w");
//				if (fp != 0) {
//					for (uint i = 0; i < msgLog().size(); ++i) {
//						fprintf(fp, "%s\n", msgLog()[i].c_str());
//					}
//					msgLog().clear();
//				} else {
//					LogToConsole(2, "LOG", "Unable to open logfile, logging to file disabled.", 0);
//					logToFile = false;
//				}
//
//				preInit = false;
//			}
//		}
	} catch (...) {
		NGMUnlock(&__Log::mutex);
		init = false;
		throw;
	}

	NGMUnlock(&__Log::mutex);
	init = false;
}

_Log const & _Log::Instance() {
	NGMOnce(&__Log::once_control, __Log::Init);

	return *__Log::pInstance;
}

void Fatal() {
	//Terminating = true;
	Log.Error("This error is fatal. Quitting...");
	ResetConsole();
	exit(1);
}

void _Log::_Message(int lvl, char const * const title, char const * const s, ...) const {
	va_list args;

	if (init) {
		printf("Log Init active - Message blocked.\n");
		printf("(lvl = %i) (t)%s %s\n", lvl, title, s);
		return;
	}
	NGMLock(&__Log::mutex);
	va_start(args, s);
	LogToConsole(lvl, title, s, args);
	va_end(args);

	/*// File Logging disabled (use "tee")
	 va_start(args, s);
	 LogToFile(lvl, title, s, args);
	 va_end(args);
	 */
	NGMUnlock(&__Log::mutex);

        if(lvl == 1) {
                warningCount += 1;
                if(warningCount > 100) {
                        printf("Max number of warnings reached!\nPlease report this issue on http://github.com/Cibiv/NextGenMap/issues!\n");
                        Fatal();
                }
        }

}

void _Log::Cleanup() {
	delete pInstance;
}

_Log::_Log() {
}
_Log::~_Log() {
}

void _Log::FilterLevel(int const lvl) {
	__Log::filterlvl = lvl;
}

