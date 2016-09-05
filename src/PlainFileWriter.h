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

#ifndef PLAINFILEWRITER_H_
#define PLAINFILEWRITER_H_

#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <cstring>

#include "NGMThreads.h"
#include "ILog.h"
#include "IConfig.h"
#include "FileWriter.h"

class PlainFileWriter: public FileWriter {

public:

	FILE * m_Output;

	PlainFileWriter(char const * const filename) {
		if (filename == 0) {
			m_Output = stdout;
		} else {
			if (!(m_Output = fopen(filename, "w"))) {
				Log.Error("Unable to open output file %s", filename);
			}
		}
	}

	~PlainFileWriter() {
		fclose(m_Output);
	}

	void doFlush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {

		if (bufferPosition > BUFFER_LIMIT || last) {
#ifdef __APPLE__
			fwrite(writeBuffer, sizeof(char), bufferPosition, m_Output);
#else
			fwrite_unlocked(writeBuffer, sizeof(char), bufferPosition, m_Output);
#endif

			bufferPosition = 0;

		}

	}

};

#endif /* PLAINFILEWRITER_H_ */
