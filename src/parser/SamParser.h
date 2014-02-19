/*
 * SamParser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include "IParser.h"

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

int const buffer_size = 10000;

class SamParser: public IParser {

private:
	gzFile fp;
	bool parse_all;
	char *buffer;

	kseq_t * tmp;

public:
	virtual ~SamParser() {
		delete[] buffer;
		buffer = 0;
		gzclose(fp);
	}

	virtual void init(char const * fileName);
	virtual size_t doParseRead(MappedRead * read);
};

#endif /* SAMPARSER_H_ */
