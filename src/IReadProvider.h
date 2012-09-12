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

	virtual uint init(char const * fileName) = 0;

	virtual MappedRead * GenerateRead(int const readid) = 0;

	virtual void DisposeRead(MappedRead * read) = 0;
};

#endif /* IREADPROVIDER_H_ */
