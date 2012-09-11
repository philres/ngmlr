#ifdef _WIN32

#include <windows.h>

#include <vector>
#include <map>

#include "PlatformSpecifics.h"
#include "Log.h"
#include "NGM.h"

#include "Debug.h"

#undef module_name
#define module_name "P:WIN"

namespace _PlatformSpecifics_windows
{
	int colors[] = {-1, 15, 6, 14, 4, 12};

	HANDLE hConsole = INVALID_HANDLE_VALUE;
	COORD origin = {0, 0};
	int defaultColor = 7;

	void _InitConsole()
	{
		if (hConsole == INVALID_HANDLE_VALUE)
		{
			hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

			CONSOLE_SCREEN_BUFFER_INFO info;
			GetConsoleScreenBufferInfo(hConsole, &info);
			defaultColor = info.wAttributes;
			origin.X = info.dwCursorPosition.X;
			origin.Y = info.dwCursorPosition.Y;
		}
	}

	struct MappedFile
	{
		HANDLE hFile;
		HANDLE hMapping;
		char const * pBase;
	};

	std::vector<MappedFile> * mappedFiles;
	std::map<int, HMODULE> * loadedDLLs;
}

using namespace _PlatformSpecifics_windows;

void InitConsole()
{
	_InitConsole();
}

void SetConsoleColor(ConsoleColor color)
{
	int n = colors[color];
	if (color == -1 || n == -1)
		ResetConsoleColor();
	else
		SetConsoleTextAttribute(hConsole, n);
}

void ResetConsoleColor()
{
	SetConsoleTextAttribute(hConsole, defaultColor);
}

void ClearConsole()
{
	SetConsoleCursorPosition(hConsole, origin);
}

void Sleep(int msec)
{
	SleepEx(msec, false);
}

bool const FileExists(char const * const filename)
{
	if (filename == 0)
	{
		Log.Warning("Invalid filename passed to FileExists");
		return false;
	}
	FILE * fp;
	fopen_s(&fp, filename, "r");
	if ( fp )
	{
		fclose(fp);
		return true;
	}
	return false;
}

int const CreateMapping(char const * const filename, char const * &pData)
{
	MappedFile mf;
	mf.hFile = CreateFile(filename, GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	mf.hMapping = CreateFileMapping(mf.hFile, 0, PAGE_READONLY, 0, 0, 0);

	if (mf.hMapping == NULL || mf.hMapping == INVALID_HANDLE_VALUE)
	{
		Log.Error("Could not create file mapping object (%d).\n", GetLastError());
		Fatal();
	}

	mf.pBase = pData = (char const *)MapViewOfFile(mf.hMapping, FILE_MAP_READ, 0, 0, 0);

	mappedFiles->push_back(mf);
	return mappedFiles->size()-1;
}

void CloseMapping(int const mapping)
{
	if (mapping >= mappedFiles->size())
	{
		Log.Warning("Tried to close unknown mapped file %i", mapping);
		return;
	}

	MappedFile mf = (*mappedFiles)[mapping];

	if (mf.pBase == 0)
	{
		Log.Warning("Tried to close already closed mapped file %i", mapping);
		return;
	}
	UnmapViewOfFile(mf.pBase);
	mf.pBase = 0;
	CloseHandle(mf.hMapping);
	CloseHandle(mf.hFile);
}

#include "ILog.h"
#include "IConfig.h"

int const InitDLL(char const * const filename)
{
	try
	{
		int err = 0;

		HMODULE sdll = LoadLibrary(filename);
		err = GetLastError();
		if (err != 0) 
			throw err;

		pfSetLog SetLog = (pfSetLog)GetProcAddress(sdll, "SetLog");
		if ((err = GetLastError()) == 0)
			SetLog(&Log);
#ifdef _DEBUG
		else
			Log.Warning("No SetLog method exported from %s", filename);
#endif

		pfSetConfig SetConfig = (pfSetConfig)GetProcAddress(sdll, "SetConfig");
		if ((err = GetLastError()) == 0)
			SetConfig(&Config);
#ifdef _DEBUG
		else
			Log.Warning("No SetConfig method exported from %s", filename);
#endif
		int pos = loadedDLLs->size();
		(*loadedDLLs)[pos] = sdll;
		return pos;
	}
	catch (int err)
	{
		Log.Error("Unable to load DLL %s (Error %i)", filename, err);
	}
	return 0;
}

void * GetDLLFunc(int const dll, char const * const name, bool required)
{
	HMODULE sdll = (*loadedDLLs)[dll];
	void * pf = GetProcAddress(sdll, name);
	if (required && (pf == 0))
	{
		Log.Error("Unable to load function %s", name);
		Fatal();
	}
	return pf;
}

void InitPlatform()
{
	mappedFiles = new std::vector<MappedFile>();
	loadedDLLs = new std::map<int, HMODULE>();
}

void CleanupPlatform()
{
	delete loadedDLLs;
	delete mappedFiles;
}

void ResetConsole()
{

}

int GetPID()
{
	return GetCurrentProcessId();
}

#endif
