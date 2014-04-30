#ifndef __ILOG_H__
#define __ILOG_H__

class ILog
{
public:
	virtual void _Message(int const lvl, char const * const title, char const * const msg, ...) const = 0;
	virtual void _Debug  (int const lvl, char const * const title, char const * const msg, ...) const = 0;
	virtual ~ILog() {};
	void * null;
};

enum {
	LOG_INFO = 1,
	LOG_INPUT = 2,
	LOG_OUTPUT = 4,
	LOG_CS = 8,
	LOG_SCORES = 16,
	LOG_ALIGN = 32,
	LOG_RESULTS_SCORES = 64,
	LOG_RESULTS_ALIGN = 128,
	LOG_BUFFER_SCORES = 256,
	LOG_BUFFER_ALIGN = 512,
	LOG_SCORE_DETAILS = 1024,
	LOG_ALIGN_DETAILS = 2048,
	LOG_CS_DETAILS = 4096,
	LOG_RESULTS_CS = 8192,
	LOG_INPUT_DETAILS = 16384,
	LOG_OUTPUT_DETAILS = 32768
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

#ifdef DEBUGLOG
#define Debug(lvl, s, ...) _Debug(lvl, module_name, s, ##__VA_ARGS__)
#else
#define Debug(lvl, s, ...) null
#endif

#endif
