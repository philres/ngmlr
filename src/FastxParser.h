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

#ifndef FASTXPARSER_H_
#define FASTXPARSER_H_

#include "IParser.h"

#include <zlib.h>

class FastXParser: public IParser {

private:
	gzFile fp;

	kseq_t * tmp;

public:

	FastXParser(int const p_qryMaxLen) : IParser(p_qryMaxLen) {
		fp = 0;
		tmp = 0;
	}

	virtual ~FastXParser() {
		if(tmp != 0)
		{
			kseq_destroy(tmp);
			tmp = 0;
		}
		gzclose(fp);
	}

	virtual void init(char const * fileName) {
		fp = gzopen(fileName, "r");
		tmp = kseq_init(fp);
	}

	virtual int doParseRead(MappedRead * read) {
		int l = kseq_read(tmp);
		return copyToRead(read, tmp, l);
	}

//	virtual int doParseRead(SAMRecord * read) {
//		int l = kseq_read(tmp);
//		return copyToRead(read, tmp, l);
//	}
};

#endif /* FASTXPARSER_H_ */
