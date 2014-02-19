/*
 * Parser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef PARSER_H_
#define PARSER_H_

#include "MappedRead.h"

#include <zlib.h>
//#include <stdio.h>
#include "kseq.h"
//
KSEQ_INIT(gzFile, gzread)

class IParser {
//protected:
//	gzFile fp;
public:

	static size_t const MAX_READNAME_LENGTH = 100;

	int qryMaxLen;

//	MappedRead * read;

	virtual ~IParser() {
//		if (read != 0) {
//			kseq_destroy(read);
//			read = 0;
//		}
	}

	virtual void init(char const * fileName) = 0;
	size_t parseRead(MappedRead * pRead) {
		assert(pRead != 0);
		return doParseRead(pRead);
	}

protected:

	virtual size_t doParseRead(MappedRead * pRead) = 0;

	void copyToRead(MappedRead * read, kseq_t * kseq) {

		int nameLength = std::min(MAX_READNAME_LENGTH - 1, kseq->name.l);
		memcpy(read->name, kseq->name.s, nameLength);
		read->name[nameLength] = '\0';

		//Sequence
		memset(read->Seq, '\0', qryMaxLen);
		if (kseq->seq.l != 0) {
			read->length = std::min(kseq->seq.l,
					(size_t) qryMaxLen - 1);
			int nCount = 0;
			for (int i = 0; i < read->length; ++i) {
				char c = toupper(kseq->seq.s[i]);
				if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
					read->Seq[i] = c;
				} else {
					read->Seq[i] = 'N';
					nCount += 1;
				}

			}
		} else {
			//Log.Verbose("Empty read found (%s). Filling with Ns.", read->name);
			read->length = qryMaxLen - 2;
			memset(read->Seq, 'N', read->length);
			read->SetFlag(NGMNames::Empty);
		}

		//Log.Message("%s", read->Seq);

		//Quality
		//read->qlty = 0;
		if (kseq->qual.l > 0) {
			memcpy(read->qlty, kseq->qual.s, read->length);
			read->qlty[read->length] = '\0';
		} else {
			read->qlty[0] = '*';
			read->qlty[1] = '\0';
		}
	}

};

#endif /* PARSER_H_ */
