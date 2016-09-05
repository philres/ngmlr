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

#ifndef GZFILEWRITER_H_
#define GZFILEWRITER_H_

#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <cstring>

#include "NGMThreads.h"
#include "ILog.h"
#include "zlib.h"

class GZFileWriter : public FileWriter {

public:

	gzFile m_Output;

	GZFileWriter(char const * const filename) {
		std::string strFilename = std::string(filename);
		if(!hasSuffix(strFilename, ".gz")) {
			strFilename += ".gz";
		}
		if (!(m_Output = gzopen(strFilename.c_str(), "wb"))) {
			Log.Error("Unable to open output file %s", filename);
		}
	}

	~GZFileWriter() {
		gzclose(m_Output);
	}

	void doFlush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {

		if (bufferPosition > BUFFER_LIMIT || last) {
			if(gzwrite(m_Output, writeBuffer, bufferPosition) < 0) {
				Log.Error("Writing");
			}
			bufferPosition = 0;

		}

	}

private:
	bool hasSuffix(const std::string &str, const std::string &suffix)
	{
	    return str.size() >= suffix.size() &&
	           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}

};

#endif /* GZFILEWRITER_H_ */
