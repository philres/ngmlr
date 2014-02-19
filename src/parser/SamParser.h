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

	bool parseAdditionalInfo;

public:

//	SamParser(bool const keepTags) :
//			parseAdditionalInfo(keepTags) {
//
//	}

	virtual ~SamParser() {
		if (tmp != 0) {
			kseq_destroy(tmp);
			tmp = 0;
		}
		delete[] buffer;
		buffer = 0;
		gzclose(fp);
	}

	virtual void init(char const * fileName, bool const keepTags);
	virtual size_t doParseRead(MappedRead * read);
};

#endif /* SAMPARSER_H_ */
