/**
 * Contact: philipp.rescheneder@gmail.com
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
