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

	virtual void init(char const * fileName){
		fp = gzopen(fileName, "r");
		read = kseq_init(fp);
	}

	virtual size_t parseRead(){
		return  kseq_read(read);
	}

};


#endif /* FASTXPARSER_H_ */
