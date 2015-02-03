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
#include <string>
#include <map>

struct VcfSNP
{
	uloc pos;
	std::string ref;
	std::string alt;
};

class VcfParser {
	
private:
	std::vector<VcfSNP> snps;
	std::map<std::string, int> refmap;

	uint getRefStart(std::string ref);
	void parse_line(std::string line, uint line_num);
	void add_line(std::string chrom, std::string pos, std::string ref, std::string alt, uint line_num);

	bool isSequence(const std::string& what);

public:
	 VcfParser();
	~VcfParser();

	void open(const char* filename);

	const VcfSNP& get(uint i) const { return snps[i]; }
	uint length() const { return snps.size(); }
	bool empty() const { return snps.size() == 0; }
};


#endif /* VCFPARSER_H_ */
