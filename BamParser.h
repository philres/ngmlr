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

	virtual kseq_t * init_seq(char const * fileName);
	virtual size_t parseRead(kseq_t *& read);

};
#endif

#endif /* BAMPARSER_H_ */
