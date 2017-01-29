/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef __TIMING2_H__
#define __TIMING2_H__

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
	}
};

#endif
