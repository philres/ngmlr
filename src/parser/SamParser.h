/*
 * SamParser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include "IParser.h"

int const buffer_size = 10000;

class SamParser: public IParser {

private:
//	gzFile fp;
	bool parse_all;
	char *buffer;
public:
	virtual ~SamParser() {
		delete[] buffer;
		buffer = 0;
		gzclose(fp);
	}

	virtual kseq_t * init_seq(char const * fileName);
	virtual size_t parseRead(kseq_t *& read);
};

#endif /* SAMPARSER_H_ */
