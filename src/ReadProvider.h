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

#ifndef READPROVIDER_H_
#define READPROVIDER_H_

#include "IReadProvider.h"
#include "IRefProvider.h"
#include "IParser.h"
#include "kseq.h"

class ReadProvider: public IReadProvider {
public:

	ReadProvider();
	virtual ~ReadProvider();

	virtual uint init();

	virtual bool GenerateRead(int const readid1, MappedRead * & read1, int const readid2, MappedRead * & read2);
	virtual void DisposeRead(MappedRead * read);

private:

	size_t const readPartLength;

	size_t const bufferLength;

	size_t parsedReads;

//	MappedRead * * readBuffer;

	size_t readsInBuffer;

	IParser * parser1;

	void splitRead(MappedRead * read);
	virtual MappedRead * NextRead(IParser * parser, int const id);
	MappedRead * GenerateSingleRead(int const readid);
	IParser * DetermineParser(char const * fileName, int const qryMaxLen);
};

#endif /* READPROVIDER_H_ */
