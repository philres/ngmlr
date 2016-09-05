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

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include "IParser.h"

#include <zlib.h>
#include <stdio.h>

#include "SAMRecord.h"
#include "kseq.h"

int const buffer_size = 10000;

class SamParser: public IParser {

private:
	gzFile fp;
	bool parse_all;
	char *buffer;

	kseq_t * tmp;

	bool parseAdditionalInfo;

public:

	SamParser(int const p_qryMaxLen) : IParser(p_qryMaxLen) {
		fp = 0;
		parse_all = true;
		buffer = 0;
		tmp = 0;
		parseAdditionalInfo = false;
	}

	virtual ~SamParser() {
		if (tmp != 0) {
			kseq_destroy(tmp);
			tmp = 0;
		}
		delete[] buffer;
		buffer = 0;
		gzclose(fp);
	}

	virtual void init(char const * fileName, bool const keepTags);
	virtual int doParseRead(MappedRead * read);
	virtual int doParseRead(SAMRecord * read);
};

#endif /* SAMPARSER_H_ */
