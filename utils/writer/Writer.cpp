/*
 * NgmWriter.cpp
 *
 *  Created on: May 4, 2011
 *      Author: philipp_
 */

#include "Writer.h"

//#include "kseq.h"
//#include <zlib.h>
//#include <stdio.h>

Writer::Writer(char const * const fileName) {

	data.open((std::string(fileName)).c_str());
	//idx.open((std::string(fileName) + ".idx").c_str(), std::ios::binary);

	//data.exceptions(ofstream::failbit | ofstream::badbit);
	//idx.exceptions(ofstream::failbit | ofstream::badbit);
}

void Writer::toUpperCase(char * sequence, int sequenceLength) {
	int i = 0;
	while (sequence[i] != '\0') {
		sequence[i] = std::toupper(sequence[i]);
		i += 1;
	}
}

Writer::~Writer() {
	data.close();
}

//void NgmWriter::WriteHeader(unsigned int seqCount, unsigned int seqFieldWith,
//		unsigned int flags, unsigned int checkSum) {
//	//if (debug) {
//	//	idx << sig << " " << headerLength << " " << seqCount << " "
//	//			<< seqFieldWith << std::endl;
//	//} else {
//	idx.write((char *) &sig, 4);
//	idx.write((char *) &version, 4);
//	idx.write((char *) &headerLength, 4);
//	idx.write((char *) &entryLength, 4);
//	idx.write((char *) &seqCount, 4);
//	idx.write((char *) &seqFieldWith, 4);
//	idx.write((char *) &flags, 4);
//	idx.write((char *) &checkSum, 4);
//	//}
//}
