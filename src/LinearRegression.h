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
