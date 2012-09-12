#ifndef __ILOG_H__
#define __ILOG_H__

class ILog
{
public:
	virtual void _Message(int const lvl, char const * const title, char const * const msg, ...) const = 0;
	virtual ~ILog() {};
	void * null;
};

typedef void (*pfSetLog)(ILog const *);

extern ILog const * _log;
#define Log (*_log)

#undef module_name
#define module_name 0

#define Message(s, ...) _Message(0, module_name, s , ##__VA_ARGS__)
#define Warning(s, ...) _Message(1, module_name, s , ##__VA_ARGS__)
#define Error(s, ...) _Message(2, module_name, s , ##__VA_ARGS__)
#define Green(s, ...) _Message(3, module_name, s , ##__VA_ARGS__)
#define Progress(s, ...) _Message(99, "Progress", s , ##__VA_ARGS__)

//#define VERBOSE

#ifdef VERBOSE
#define Verbose(s, ...) _Message(0, module_name, s, ##__VA_ARGS__)
#else
#define Verbose(s, ...) null
#endif

#endif
