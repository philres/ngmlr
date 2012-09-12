/*
 * VcfParser.h
 *
 *  Created on: Sep 11, 2012
 *      Author: fritz
 */

#ifndef VCFPARSER_H_
#define VCFPARSER_H_
#include <map>

struct SNP{
	unsigned int pos;
	char alt;
};

int const buffer_size = 300;

class VcfParser{

private:
	char *buffer;
	gzFile fp;
	std::list <SNP> snp_list;
	bool VcfParser::compare_SNP(SNP first, SNP second);

public:
	~VcfParser(){
		gzclose(fp);
		delete buffer;
		snp_list.clear();
	}

	VcfParser(char const * fileName);
	SNP getNextSnp();
};


#endif /* VCFPARSER_H_ */
