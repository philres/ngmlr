/*
 * NgmWriter.h
 *
 *  Created on: May 4, 2011
 *      Author: philipp_
 */

#ifndef NGMWRITER_H_
#define NGMWRITER_H_

#include <fstream>
#include <string.h>

#include "MappedRead.h"

using std::ofstream;

//struct kseq_t;

//KSEQ_INIT(gzFile, gzread)

class Writer {

protected:
	ofstream data;

//	void WriteHeader(unsigned int seqCount,
//			unsigned int seqFieldWith, unsigned int flags,
//			unsigned int checkSum);

	void toUpperCase(char * sequence, int sequenceLength);

public:
	Writer(char const * const fileName);

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
