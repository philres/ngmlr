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

	IParser * parser1;
	IParser * parser2;

	char const peDelimiter;

	bool const isPaired;

	bool const skipMateCheck;

	virtual MappedRead * NextRead(IParser * parser, int const id);
	MappedRead * GenerateSingleRead(int const readid);
	IParser * DetermineParser(char const * fileName, int const qryMaxLen);
};

#endif /* READPROVIDER_H_ */
