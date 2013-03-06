/*
 * FileWriter.h
 *
 *  Created on: Mar 2, 2013
 *      Author: philipp
 */

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include <stdio.h>
#include <stdarg.h>

#include "NGMThreads.h"
#include "ILog.h"
#include "Config.h"

#include <iostream>
#include <cstring>

static int count = 0;

class FileWriter {

public:

	FILE * m_Output;

	NGMMutex m_OutputMutex;

	FileWriter(char const * const filename) {
		count += 1;
		if(count > 1) {
			Log.Error("To many FileWriter instances!");
			Fatal();
		}
		NGMInitMutex(&m_OutputMutex);
		if (!(m_Output = fopen(filename, "w"))) {
			Log.Error("Unable to open output file %s", filename);
			Fatal();
		}

	}

	~FileWriter() {
		fclose(m_Output);
	}

	void Flush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {
		//NGM.AquireOutputLock();
		//		Log.Error("Lock");
		NGMLock(&m_OutputMutex);
		if (bufferPosition > BUFFER_LIMIT || last) {
//			NGM.AquireOutputLock();

			//std::cerr << std::string(writeBuffer, bufferPosition) << std::endl;
			//Log.Message("%.*s, bufferPosition, writeBuffer");
			if (fwrite_unlocked(writeBuffer, sizeof(char), bufferPosition, m_Output)
					< 0) {
				Log.Error("Writing");
				Fatal();
			}
			fflush(m_Output);
//			NGM.ReleaseOutputLock();
			bufferPosition = 0;

		}
		NGMUnlock(&m_OutputMutex);
		//		Log.Green("Unlock");
		//NGM.ReleaseOutputLock();

	}

};

#endif /* FILEWRITER_H_ */
