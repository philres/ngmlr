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

#ifndef IREADPROVIDER_H_
#define IREADPROVIDER_H_

#include "MappedRead.h"

class IReadProvider {
public:

	virtual ~IReadProvider() { }

	virtual uint init() = 0;

	virtual bool GenerateRead(int const readid1, MappedRead * & read1, int const readid2, MappedRead * & read2) = 0;

	virtual void DisposeRead(MappedRead * read) = 0;
};

#endif /* IREADPROVIDER_H_ */
