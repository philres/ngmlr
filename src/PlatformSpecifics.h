#ifndef __PLATFORMSPECIFICS_H__
#define __PLATFORMSPECIFICS_H__

#include "Types.h"

enum ConsoleColor
{
	Message,
	MessageTitle,
	Warning,
	WarningTitle,
	Error,
	ErrorTitle
};

void InitPlatform();
void CleanupPlatform();

void InitConsole();
void ClearConsole();

void SetConsoleColor(ConsoleColor const color);
void ResetConsoleColor();
void ResetConsole();

#ifndef _WIN32
ulong GetTickCount();
#endif

int GetPID();

// Blocking wait for msec milliseconds
void Sleep(int msec);

bool const FileExists(char const * const filename);
uloc const FileSize(char const * const filename);
int const CreateMapping(char const * const filename, char const * &pData);
void Remap(int const mapping, char const * & pData);
void CloseMapping(int const mapping);
uloc GetMapLength(int const map);

//bool CheckDLL(char const * const filename);
//int const InitDLL(char const * const filename);
//void * GetDLLFunc(int const dll, char const * const name, bool required = true);

#ifndef _WIN32
int _kbhit(void);
int _getch();
#endif

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

// Atomically increment *pi and return the incremented value
// _asm: LOCK INC DWORD PTR pi
inline int AtomicInc(volatile int * pi)
{
#ifdef _WIN32
	return (int)InterlockedIncrement((volatile long*)pi);
#endif

#ifndef _WIN32
	// GCC Builtin
	return __sync_add_and_fetch(pi, 1);
#endif
}

inline int AtomicDec(volatile int * pi)
{
#ifdef _WIN32
	return (int)InterlockedDecrement((volatile long*)pi);
#endif

#ifndef _WIN32
	return __sync_add_and_fetch(pi, -1);
#endif
}

#endif
