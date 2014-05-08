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

class GZFileWriter : public FileWriter {

public:

	gzFile m_Output;

	GZFileWriter(char const * const filename) {
		std::string strFilename = std::string(filename) + ".gz";
		if (!(m_Output = gzopen(strFilename.c_str(), "wb"))) {
			Log.Error("Unable to open output file %s", filename);
			Fatal();
		}
	}

	~GZFileWriter() {
		gzclose(m_Output);
	}

	void doFlush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {

		if (bufferPosition > BUFFER_LIMIT || last) {
			if(gzwrite(m_Output, writeBuffer, bufferPosition) < 0) {
				Log.Error("Writing");
				Fatal();
			}
			bufferPosition = 0;

		}

	}

};

#endif /* GZFILEWRITER_H_ */
