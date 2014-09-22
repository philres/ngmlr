/*
 * FileWriter.h
 *
 *  Created on: Mar 2, 2013
 *      Author: philipp
 */

#ifndef PLAINFILEWRITER_H_
#define PLAINFILEWRITER_H_

#include <stdio.h>
#include <stdarg.h>

#include "NGMThreads.h"
#include "ILog.h"
#include "Config.h"
#include "FileWriter.h"

#include <iostream>
#include <cstring>


class PlainFileWriter : public FileWriter {

public:

	FILE * m_Output;

	PlainFileWriter(char const * const filename) {
		if (!(m_Output = fopen(filename, "w"))) {
			Log.Error("Unable to open output file %s", filename);
			Fatal();
		}
	}

	~PlainFileWriter() {
		fclose(m_Output);
	}

	void doFlush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {

		if (bufferPosition > BUFFER_LIMIT || last) {
			if (fwrite_unlocked(writeBuffer, sizeof(char), bufferPosition, m_Output) < 0) {
				Log.Error("Writing");
				Fatal();
			}
			bufferPosition = 0;

		}

	}

};

#endif /* PLAINFILEWRITER_H_ */
