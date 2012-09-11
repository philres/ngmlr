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


class FastXParser:public IParser{

private:
	gzFile fp;

public:

	virtual ~FastXParser() {
		gzclose(fp);
	}

	virtual kseq_t * init_seq(char const * fileName){
		fp = gzopen(fileName, "r");
		kseq_t *seq = kseq_init(fp);
		return seq;
	}

	virtual size_t parseRead(kseq_t *& read){
		return  kseq_read(read);
	}

};


#endif /* FASTXPARSER_H_ */
