#ifndef __NGMTHREADS_H__
#define __NGMTHREADS_H__

// Simple Windows/pThread wrapper

#ifdef _WIN32
// Windows Threads
#define NOMINMAX
#include <windows.h>

#define NGMTHREADFUNC DWORD WINAPI

typedef CRITICAL_SECTION NGMMutex;
typedef HANDLE NGMThread;
typedef HANDLE NGMThreadWait;

struct NGMOnceControl
{
	bool done;
	LPCRITICAL_SECTION cs;
	NGMOnceControl();
	~NGMOnceControl();
};
#define NGM_ONCE_INIT NGMOnceControl();

#endif
#ifndef _WIN32
// POSIX Threads
#include "pthread.h"
#include "sched.h"

#define NGMTHREADFUNC void*

typedef pthread_t NGMThread;
typedef pthread_cond_t NGMThreadWait;
typedef pthread_mutex_t NGMMutex;

typedef pthread_once_t NGMOnceControl;
#define NGM_ONCE_INIT PTHREAD_ONCE_INIT

#endif

typedef NGMTHREADFUNC NGMThreadFunc(void*);
typedef NGMThreadFunc* pNGMThreadFunc;

NGMThread NGMCreateThread(pNGMThreadFunc func, void * data, bool detached = false);
void NGMJoin(NGMThread * thread);
void NGMJoinAll(NGMThread * thread, int count);

void NGMSetThreadAffinity(NGMThread * thread, int cpu);

void NGMInitWait(NGMThreadWait * wait);
void NGMWait(NGMMutex * mutex, NGMThreadWait * wait);
void NGMSignal(NGMThreadWait * wait);

void NGMInitMutex(NGMMutex * mutex);
void NGMLock(NGMMutex * mutex);
void NGMUnlock(NGMMutex * mutex);

void NGMOnce(NGMOnceControl * once_control, void (*init_func)());

#endif
