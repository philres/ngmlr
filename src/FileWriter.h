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

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include "NGMThreads.h"

class FileWriter {

public:

	NGMMutex m_OutputMutex;

	FileWriter() {
		NGMInitMutex(&m_OutputMutex);
	}

	virtual ~FileWriter() {
	}

	void Flush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {
		NGMLock(&m_OutputMutex);
		doFlush(bufferPosition, BUFFER_LIMIT, writeBuffer, last);
		NGMUnlock(&m_OutputMutex);
	}

protected:
	virtual void doFlush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) = 0;

};

#endif /* FILEWRITER_H_ */
