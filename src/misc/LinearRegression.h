#ifndef __LINEAR_REGRESSION_H__
#define __LINEAR_REGRESSION_H__

#include <stdlib.h>
#include <math.h>                           /* math functions                */


//#define REAL float
#define REAL double

inline static REAL sqr(REAL x) {
	return x * x;
}

int linreg(int n, const REAL x[], const REAL y[], REAL* m, REAL* b, REAL* r);


#endif //__LINEAR_REGRESSION_H__
