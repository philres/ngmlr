#ifdef _WIN32

#include "NGM.h"

#include "Debug.h"

//NGMThread NGMCreateThread(pNGMThreadFunc func, int tid)
//{
//	return NGMCreateThread(func, new int(tid));
//}
NGMThread NGMCreateThread(pNGMThreadFunc func, void* data, bool detached)
{
	return CreateThread(0, 0, func, data, 0, 0);
}

void NGMJoin(NGMThread * thread)
{
	WaitForSingleObject(*thread, INFINITE);
}
void NGMJoinAll(NGMThread * threads, int count)
{
	WaitForMultipleObjects(count, threads, true, INFINITE);
}

void NGMInitWait(NGMThreadWait * wait)
{
	*wait = CreateEvent(0, false, false, 0);
}
void NGMWait(NGMMutex * mutex, NGMThreadWait * wait)
{
	LeaveCriticalSection(mutex);
	WaitForSingleObject(*wait, INFINITE);
	EnterCriticalSection(mutex);
}
void NGMSignal(NGMThreadWait * wait)
{
	SetEvent(*wait);
}

void NGMInitMutex(NGMMutex * mutex)
{
	InitializeCriticalSection(mutex);
}
void NGMLock(NGMMutex * mutex)
{
	EnterCriticalSection(mutex);
}
void NGMUnlock(NGMMutex * mutex)
{
	LeaveCriticalSection(mutex);
}

NGMOnceControl::NGMOnceControl()
{
	done = false;
	cs = new CRITICAL_SECTION();
	try
	{
		InitializeCriticalSection(cs);
	}
	catch (int ex)
	{
		if (ex == STATUS_NO_MEMORY)
			Log.Error("Unable to initialize once_control: Out of memory?");
	}
}

NGMOnceControl::~NGMOnceControl()
{
	delete cs;
}

void NGMOnce(NGMOnceControl * once_control, void (*init_func)())
{
	if (!once_control->done)
	{
		EnterCriticalSection(once_control->cs);
		if (!once_control->done)
			init_func();
		once_control->done = true;
		LeaveCriticalSection(once_control->cs);
	}
}

void NGMSetThreadAffinity(NGMThread * thread, int cpu)
{
	HANDLE hSelf = GetCurrentThread();
	if (thread == 0)
		thread = &hSelf;
	SetThreadAffinityMask(*thread, 1 << cpu);
}

#endif
