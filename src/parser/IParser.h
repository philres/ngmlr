/*
 * Parser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef PARSER_H_
#define PARSER_H_
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

class IParser {
protected:
	gzFile fp;
public:

	virtual ~IParser() {
		if (read != 0) {
			kseq_destroy(read);
			read = 0;
		}
	}
	virtual void init(char const * fileName) = 0;
	virtual size_t parseRead() = 0;

	kseq_t * read;
};

#endif /* PARSER_H_ */
