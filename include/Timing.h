#ifndef __TIMING_H__
#define __TIMING_H__

#include <time.h>

#ifndef _WIN32
#include <sys/time.h>

class Timer
{
private:
	timeval start, end;

public:
	inline void ST()
	{
		gettimeofday(&start, 0);
	}
	inline volatile float ET()
	{
		gettimeofday(&end, 0);

		return (float)(end.tv_sec - start.tv_sec) + (float)(end.tv_usec - start.tv_usec) / 1000000.0f;
	}
};

#endif
#ifdef _WIN32

class Timer
{
private:
    clock_t clock_start;

public:
	inline void ST()
	{
		clock_start = clock();
	}
	inline volatile float ET()
	{
		clock_t diff = clock();
		diff -= clock_start;
		return (float)diff/(float)CLOCKS_PER_SEC;
	}
};

#endif


#endif
