#ifndef _WIN32

#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <dlfcn.h>
#include "Log.h"

#include <map>
#include <termios.h>
#include <memory.h>

#include <time.h>

#include <signal.h>

#include <stdlib.h>

#include "PlatformSpecifics.h"

#include <errno.h>

#undef module_name
#define module_name "P:UNIX"

namespace _PlatformSpecifics_unix
{
	char const * const colors[] = { "0", "1;30", "0;33", "1;33", "0;31", "1;31", "0;32", "1;32" };

	std::map<int, void *> & loadedDLLs()
	{
		static std::map<int, void *> * res = new std::map<int, void *>();
		return *res;
	}

	std::map<int, std::pair<char const *, uloc> > & mappings()
	{
		static std::map<int, std::pair<char const *, uloc> > * res = new std::map<int, std::pair<char const *, uloc> >();
		return *res;
	}
}

using namespace _PlatformSpecifics_unix;

void SetConsoleColor(ConsoleColor color) {
	fprintf(stderr, "\e[%sm", colors[color]);
}

void ResetConsoleColor() {
	fprintf(stderr, "\e[0m");
}

void ClearConsole() {
}

void Sleep(int msec) {
	timespec tme;

	int sec = 0;
	if (msec > 1000)
	{
		sec = msec / 1000;
		msec -= sec*1000;
	}
	tme.tv_sec = sec;
	tme.tv_nsec = msec * 1000000;
	nanosleep(&tme, 0);
}

int GetPID()
{
	return getpid();
}

bool const FileExists(char const * const filename) {
	FILE * fp = fopen(filename, "r");
	if (fp) {
		fclose(fp);
		return true;
	}
	return false;
}

int const CreateMapping(char const * const filename, char const * &pData) {
	int fd = open(filename, O_RDONLY);
	if (fd == -1)
	{
		Log.Error("Error opening file %s for mapping (error %i).", filename, errno);
		Fatal();
	}
	uloc len = FileSize(filename);

#ifdef __APPLE__
	pData = (const char*) mmap(0, len, PROT_READ, MAP_PRIVATE, fd, 0);
#else
	pData = (const char*) mmap64(0, len, PROT_READ, MAP_PRIVATE, fd, 0);
#endif
	if (pData == MAP_FAILED)
	{
		Log.Error("Error mapping file %s into memory (Error returned from mmap: %i)", filename, errno);
		if (errno == EINVAL)
		{
			Log.Error("Args: Addr=0, Len=%llu, Offset=0", len);
		}
		else if (errno == ENOMEM)
		{
			Log.Error("Try \"ulimit -v\"");
		}
		Fatal();
	}

	mappings()[fd] = std::pair<char const*, uloc>(pData, len);

	return fd;
}
uloc GetMapLength(int const map)
{
	if (mappings().count(map) != 1)
	{
		Log.Error("Couldnt get size of invalid mapping %i", map);
		return 0;
	}
	else
		return mappings()[map].second;

}
void Remap(int const mapping, char const * & pData)
{
	if (mappings().count(mapping) != 1)
	{
		Log.Error("Remap error: Invalid mapping %i", mapping);
		return;
	}

	void const * addr = mappings()[mapping].first;
	uloc len = mappings()[mapping].second;

	int ret = munmap((void*) addr, len); //TODO_GENOMESIZE: munmap64?

	if (ret != 0)
		Log.Error("Remap error: Unmap returned %i (addr = 0x%x, len = %llu)", errno, addr, len);
#ifdef __APPLE__
	pData = (const char*)mmap(0, len, PROT_READ, MAP_PRIVATE, mapping, 0);
#else
	pData = (const char*)mmap64(0, len, PROT_READ, MAP_PRIVATE, mapping, 0);
#endif
	mappings()[mapping].first = pData;

	if (pData == MAP_FAILED)
	{
		Log.Error("Error remapping %i (Error returned from mmap: %i)", mapping, errno);
		Fatal();
	}
}
void CloseMapping(int const mapping) {
	munmap((void*) mappings()[mapping].first, mappings()[mapping].second);
	close(mapping);
}

//bool CheckDLL(char const * const filename) {
//	Log.Verbose("Checking dll: %s", filename);
//	return dlopen(filename, RTLD_LAZY) != 0;
//}
//
//int const InitDLL(char const * const filename) {
//	try {
//		void * sdll = dlopen(filename, RTLD_LAZY);
//
//		if (!sdll) {
//			Log.Error("Cannot open library: %s", dlerror());
//			return -1;
//		}
//
//		int pos = loadedDLLs().size();
//		loadedDLLs()[pos] = sdll;
//		pfSetLog SetLog = (pfSetLog) GetDLLFunc(pos, "SetLog");
//		if (SetLog != 0) {
//			SetLog(&Log);
//		}
//#ifdef _DEBUG
//		else
//		Log.Warning("No SetLog method exported from %s", filename);
//#endif
//
//		pfSetConfig SetConfig = (pfSetConfig) GetDLLFunc(pos, "SetConfig");
//		if (SetConfig != 0) {
//			SetConfig(&Config);
//		}
//#ifdef _DEBUG
//		else
//		Log.Warning("No SetConfig method exported from %s", filename);
//#endif
//
//		return pos;
//	} catch (int err) {
//		Log.Error("Unable to load DLL %s (Error %i)", filename, err);
//	}
//	return -1;
//}

//void * GetDLLFunc(int const dll, char const * const name, bool required) {
//	void * sdll = loadedDLLs()[dll];
//	void * pf = dlsym(sdll, name);
//	if (required && (pf == 0))
//	{
//		Log.Error("Unable to load function %s", name);
//		Fatal();
//	}
//	return pf;
//}

struct termios alt;

int _getch() {
	static int ch = -1, fd = 0;
	struct termios neu;
	fd = fileno(stdin);
	tcgetattr(fd, &alt);
	neu = alt;
	neu.c_lflag &= ~(ICANON | ECHO);
	tcsetattr(fd, TCSANOW, &neu);
	ch = getchar();
	tcsetattr(fd, TCSANOW, &alt);
	return ch;
}

int _kbhit(void) {
	struct termios term, oterm;
	int fd = 0;
	int c = 0;
	tcgetattr(fd, &oterm);
	memcpy(&term, &oterm, sizeof(term));
	term.c_lflag = term.c_lflag & (!ICANON);
	term.c_cc[VMIN] = 0;
	term.c_cc[VTIME] = 1;
	tcsetattr(fd, TCSANOW, &term);
	c = getchar();
	tcsetattr(fd, TCSANOW, &oterm);
	if (c != -1)
		ungetc(c, stdin);
	return ((c != -1) ? 1 : 0);
}

void InitConsole()
{
	tcgetattr(STDIN_FILENO, &alt);
}

void ResetConsole()
{
	tcsetattr(STDIN_FILENO, TCSANOW, &alt);
}

void HandleSignal(int sig)
{
	Log.Error("System sent signal %i", sig);
	Fatal();
}

void MuteSignal(int sig)
{
	ResetConsoleColor();
	fprintf(stderr, "Resources freed by OS. Exiting normally...\n");
	exit(0);
}

void InitPlatform()
{
	/*void (*prev_fn)(int);
	prev_fn = signal(SIGSEGV, HandleSignal);
	Log.Warning("SIGSEGV Intercept enabled.");*/
}

void CleanupPlatform()
{
	signal(SIGSEGV, MuteSignal);
}

ulong GetCurrentTicks()
{
	timeval tv;
	gettimeofday(&tv, 0);

	return (tv.tv_sec * 1000ul + tv.tv_usec / 1000ul);
}

const ulong programStart = GetCurrentTicks();

ulong GetTickCount()
{
	return GetCurrentTicks() - programStart;
}


#endif
