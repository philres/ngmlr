/**
 * Contact: philipp.rescheneder@gmail.com
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
