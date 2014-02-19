/*
 * BamParser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef BAMPARSER_H_
#define BAMPARSER_H_

#ifdef _BAM

#include "IParser.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

using namespace BamTools;

class BamParser: public IParser {

private:
	gzFile fp;
	BamMultiReader reader;
	bool parse_all;
	BamAlignment* al;
	kseq_t * tmp;
public:
	virtual ~BamParser() {
		reader.Close();
		delete al;
	}

	virtual void init(char const * fileName);
	virtual size_t doParseRead(MappedRead * read);

};
#endif

#endif /* BAMPARSER_H_ */
