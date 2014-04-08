/*
 * GZFileWriter.h
 *
 *  Created on: Apr 2, 2014
 *      Author: philipp_
 */

#ifndef GZFILEWRITER_H_
#define GZFILEWRITER_H_

#include <stdio.h>
#include <stdarg.h>

#include "NGMThreads.h"
#include "ILog.h"
#include "Config.h"

#include <iostream>
#include <cstring>
#include "zlib.h"

static int count = 0;

class FileWriter {

public:

	gzFile m_Output;

	NGMMutex m_OutputMutex;

	FileWriter(char const * const filename) {
		count += 1;
		if(count > 1) {
			Log.Error("To many FileWriter instances!");
			Fatal();
		}
		NGMInitMutex(&m_OutputMutex);
		if (!(m_Output = gzopen(filename, "wb"))) {
			Log.Error("Unable to open output file %s", filename);
			Fatal();
		}

	}

	~FileWriter() {
		gzclose(m_Output);
	}

	void Flush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {

		NGMLock(&m_OutputMutex);
		if (bufferPosition > BUFFER_LIMIT || last) {
			if(gzwrite(m_Output, writeBuffer, bufferPosition) < 0) {
				Log.Error("Writing");
				Fatal();
			}
			bufferPosition = 0;

		}
		NGMUnlock(&m_OutputMutex);
	}

};

#endif /* GZFILEWRITER_H_ */
