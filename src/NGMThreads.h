/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
 */

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
