/**
 * Contact: philipp.rescheneder@gmail.com
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
