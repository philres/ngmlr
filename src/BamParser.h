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

#ifndef BAMPARSER_H_
#define BAMPARSER_H_

#include "IParser.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

using namespace BamTools;

class BamParser: public IParser {

private:
	gzFile fp;
	BamMultiReader reader;

	BamAlignment* al;
	kseq_t * tmp;

	char additionalInfo[100000];

	std::vector<BamRegion> regions;

	std::map<std::string, int> readNames;

	int parsedReads;

public:

	BamParser(int const p_qryMaxlen) : IParser(p_qryMaxlen) {
		fp = 0;
		al = 0;
		tmp = 0;
		parsedReads = 0;
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

	void parseSnifflesFile(char const * fileName);
	void parseBedFile(char const * fileName);
	void parseRealignFile(char const * fileName);

	int doParseSingleRead(MappedRead * read, BamAlignment * al);

	virtual void init(char const * fileName);
	virtual int doParseRead(MappedRead * read);
	virtual int doParseRead(SAMRecord * read);
};

#endif /* BAMPARSER_H_ */
