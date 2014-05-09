/*
 * FileWriter.h
 *
 *  Created on: Mar 2, 2013
 *      Author: philipp
 */

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include "NGMThreads.h"

//static int count = 0;

class FileWriter {

public:

	NGMMutex m_OutputMutex;

	FileWriter() {
//		count += 1;
//		if (count > 1) {
//			Log.Error("To many FileWriter instances!");
//			Fatal();
//		}
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
