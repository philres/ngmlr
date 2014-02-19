/*
 * FastxParser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef FASTXPARSER_H_
#define FASTXPARSER_H_

#include "IParser.h"
#include <zlib.h>

class FastXParser: public IParser {

private:
	gzFile fp;

	kseq_t * tmp;

public:

	virtual ~FastXParser() {
		gzclose(fp);
	}

	virtual void init(char const * fileName) {
		fp = gzopen(fileName, "r");
		tmp = kseq_init(fp);
	}

	virtual size_t doParseRead(MappedRead * read) {
		int l = kseq_read(tmp);
		copyToRead(read, tmp);
		return l;
	}

};

#endif /* FASTXPARSER_H_ */
