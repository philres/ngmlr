/*
 * NgmWriter.h
 *
 *  Created on: May 4, 2011
 *      Author: philipp_
 */

#ifndef NGMWRITER_H_
#define NGMWRITER_H_

#include <stdio.h>

#include "MappedRead.h"


//struct kseq_t;

//KSEQ_INIT(gzFile, gzread)

class Writer {

private:
	FILE * m_Output;

	static int const BUFFER_SIZE = 17000000;
	static int const BUFFER_LIMIT = 16000000;

	char * writeBuffer;
	int bufferPosition;

protected:

//	void WriteHeader(unsigned int seqCount,
//			unsigned int seqFieldWith, unsigned int flags,
//			unsigned int checkSum);

	void toUpperCase(char * sequence, int sequenceLength);

public:
	Writer(char const * const fileName);

	int Print(const char *format, ...);

	void Flush(bool last = false);

	virtual void writeRead(MappedRead * read) = 0;

	void writeDummy(int const length) {
//		char * dummy = new char[length];
//		Entry read;
//		read.sequence = dummy;
//		memset(dummy, 'N', length);
//		writeRead(&read, length);
//		delete[] dummy;
//		dummy = 0;
	}

	virtual ~Writer();
};

#endif /* NGMWRITER_H_ */
