#ifndef __TIMING2_H__
#define __TIMING2_H__


//#ifndef _WIN32
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>

class Timer
{
private:
	timeval start, end;

public:
	void ST()
	{
		gettimeofday(&start, NULL);
	}
	inline volatile float ET()
	{
		gettimeofday(&end, NULL);

		return (float)(end.tv_sec - start.tv_sec) + (float)(end.tv_usec - start.tv_usec) / 1000000.0f;
//		return 0.0f;
	}
};

//#endif
//#ifdef _WIN32
//
//class Timer
//{
//private:
//    clock_t clock_start;
//
//public:
//	inline void ST()
//	{
//		clock_start = clock();
//	}
//	inline volatile float ET()
//	{
//		clock_t diff = clock();
//		diff -= clock_start;
//		return (float)diff/(float)CLOCKS_PER_SEC;
//	}
//};
//
//#endif
//
//
#endif
