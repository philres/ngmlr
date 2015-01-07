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
	std::vector<SNP> snps;

public:
	 VcfParser();
	~VcfParser();

	void open(const char* filename);

	const SNP& get( uint i ) const { return snps[i]; }
	uint length() const { return snps.size(); }
	bool empty() const { return snps.size() == 0; }
};


#endif /* VCFPARSER_H_ */
