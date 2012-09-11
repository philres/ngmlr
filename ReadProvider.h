/*
 * ReadProvider.h
 *
 *  Created on: Jun 14, 2012
 *      Author: philipp_
 */

#ifndef READPROVIDER_H_
#define READPROVIDER_H_

#include "IReadProvider.h"
#include "IRefProvider.h"

class ReadProvider: public IReadProvider {
public:

	ReadProvider();
	virtual ~ReadProvider();

	virtual uint init(char const * fileName);

	virtual MappedRead * GenerateRead(int const readid);
	virtual void DisposeRead(MappedRead * read);

private:

	virtual MappedRead * NextRead(int const id);
	MappedRead * GenerateSingleRead(int const readid);
	void DetermineParser(char const * fileName);
};

#endif /* READPROVIDER_H_ */
