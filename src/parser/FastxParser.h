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
		if(tmp != 0)
		{
			kseq_destroy(tmp);
			tmp = 0;
		}
		gzclose(fp);
	}

	virtual void init(char const * fileName, bool const keepTags) {
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
