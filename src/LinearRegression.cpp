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

#include "LinearRegression.h"

#include <stdlib.h>
#include <math.h>                           /* math functions                */

int linreg(int n, const REAL x[], const REAL y[], REAL* m, REAL* b, REAL* r) {
	REAL sumx = 0.0; /* sum of x                      */
	REAL sumx2 = 0.0; /* sum of x**2                   */
	REAL sumxy = 0.0; /* sum of x * y                  */
	REAL sumy = 0.0; /* sum of y                      */
	REAL sumy2 = 0.0; /* sum of y**2                   */

	for (int i = 0; i < n; i++) {
		sumx += x[i];
		sumx2 += sqr(x[i]);
		sumxy += x[i] * y[i];
		sumy += y[i];
		sumy2 += sqr(y[i]);
	}

	REAL denom = (n * sumx2 - sqr(sumx));
	if (denom == 0) {
		// singular matrix. can't solve the problem.
		*m = 0;
		*b = 0;
		if (r)
			*r = 0;
		return 1;
	}

	*m = (n * sumxy - sumx * sumy) / denom;
	*b = (sumy * sumx2 - sumx * sumxy) / denom;
	if (r != NULL) {
		*r = (sumxy - sumx * sumy / n) / /* compute correlation coeff     */
		sqrt((sumx2 - sqr(sumx) / n) * (sumy2 - sqr(sumy) / n));
	}

	return 0;
}

