/*
 * VcfParser.h
 *
 * Created originally on: Sep 11, 2012 by Fritz
 * Rewrite: December 19, 2014
 *      Author: moritz
 */

#ifndef VCFPARSER_H_
#define VCFPARSER_H_
#include "Types.h"
#include <vector>

struct SNP
{
	uloc pos;
	char alt;
};

class VcfParser {

private:
	uint next_i;
	std::vector<SNP> snps;

public:
	 VcfParser(const char* filename);
	~VcfParser();

	SNP getNextSnp();
};


#endif /* VCFPARSER_H_ */
