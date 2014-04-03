#ifndef __LOGGING_H__
#define __LOGGING_H__

#include "ILog.h"

namespace __Log
{
	void Init();
}

class _Log : public ILog
{
public:
	static _Log const & Instance();
	static void Init(char const * logFile, int logLvl);
	static void FilterLevel(int const lvl);
	static void setColor(bool const color);
	static void Cleanup();

	void _Message(int const lvl, char const * const title, char const * const msg, ...) const;
	void _Debug  (int const lvl, char const * const title, char const * const msg, ...) const;

private:

	_Log();
	~_Log();

	friend void __Log::Init();
};

void Fatal();

#undef Log
#define Log _Log::Instance()
#define Logger _Log::Instance()

#endif
