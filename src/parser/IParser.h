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

	virtual ~IParser() { }
	virtual kseq_t * init_seq(char const * fileName) = 0;
	virtual size_t parseRead(kseq_t *& read) = 0;
};


#endif /* PARSER_H_ */
