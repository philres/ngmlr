/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
 */

#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include <vector>
#include <string>
#include <map>

#include "Types.h"

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
