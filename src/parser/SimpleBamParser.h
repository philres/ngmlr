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

	bool parseAdditionalInfo;

	char additionalInfo[100000];

public:

	BamParser(int const p_qryMaxlen) : IParser(p_qryMaxlen) {
		fp = 0;
		parse_all = true;
		al = 0;
		tmp = 0;
		parseAdditionalInfo = false;
	}

	virtual ~BamParser() {
		if (tmp != 0) {
			kseq_destroy(tmp);
			tmp = 0;
		}
		reader.Close();
		gzclose(fp);

		delete al;
	}

	virtual void init(char const * fileName, bool const keepTags);
	virtual int doParseRead(MappedRead * read);
	virtual int doParseRead(SAMRecord * read);
};
#endif

#endif /* BAMPARSER_H_ */
