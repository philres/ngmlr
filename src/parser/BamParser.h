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
	BamMultiReader reader;
	bool parse_all;
	BamAlignment* al;
public:
	virtual ~BamParser() {
		reader.Close();
		delete al;
	}

	virtual void init(char const * fileName);
	virtual size_t parseRead();

};
#endif

#endif /* BAMPARSER_H_ */
