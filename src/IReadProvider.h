/*
 * IReadProvider.h
 *
 *  Created on: Jun 13, 2012
 *      Author: philipp_
 */

#ifndef IREADPROVIDER_H_
#define IREADPROVIDER_H_

#include "MappedRead.h"

class IReadProvider {
public:

	virtual ~IReadProvider() { }

	virtual uint init() = 0;

	virtual bool GenerateRead(int const readid1, MappedRead * & read1, int const readid2, MappedRead * & read2) = 0;

	virtual void DisposeRead(MappedRead * read) = 0;
};

#endif /* IREADPROVIDER_H_ */
