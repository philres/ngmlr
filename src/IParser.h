/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef PARSER_H_
#define PARSER_H_

#include "MappedRead.h"

#include <zlib.h>
#include <string.h>
#include <iostream>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

class IParser {

public:

	static size_t const MAX_READNAME_LENGTH = 250;

	IParser(int const qrymaxlen) :
			qryMaxLen(qrymaxlen) {

	}

	virtual ~IParser() {

	}

	virtual void init(char const * fileName) = 0;

	int parseRead(MappedRead * pRead) {
		assert(pRead != 0);
		return doParseRead(pRead);
	}
//	int parseSAMRecord(SAMRecord * pRead) {
//		assert(pRead != 0);
//		return doParseRead(pRead);
//	}

protected:

	int const qryMaxLen;

	virtual int doParseRead(MappedRead * pRead) = 0;
//	virtual int doParseRead(SAMRecord * pRead) = 0;

	int copyToRead(MappedRead * read, kseq_t * kseq, int const l) {
		int nameLength = 0;
		if (l >= 0) {
			if (kseq->seq.l == kseq->qual.l || kseq->qual.l == 0) {

				nameLength = std::min(MAX_READNAME_LENGTH - 1, kseq->name.l);
				memcpy(read->name, kseq->name.s, nameLength);
				read->name[nameLength] = '\0';

				//Sequence
				if (kseq->seq.l != 0) {
//					read->length = std::min(kseq->seq.l,
//							(size_t) qryMaxLen - 1);
					read->length = kseq->seq.l;
					read->Seq = new char[std::max(qryMaxLen, read->length + 1)];
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
					read->Seq[read->length] = '\0';
				} else {
					fprintf(stderr, "Empty read found!\n");

					read->length = qryMaxLen;
					read->Seq = new char[read->length];
					memset(read->Seq, 'N', read->length - 1);
					read->Seq[read->length] = '\0';
					read->SetFlag(NGMNames::Empty);
				}

				//Quality
				if (kseq->qual.l > 0) {
					read->qlty = new char[read->length + 1];
					memcpy(read->qlty, kseq->qual.s, read->length);
					read->qlty[read->length] = '\0';
				} else {
					read->qlty = new char[2];
					read->qlty[0] = '*';
					read->qlty[1] = '\0';
				}
			} else {
				throw "Error while parsing. Read length not equal to length of quality values!";
				//Log.Error("Discarding read %s. Length of read not equal length of quality values.", parser->read->name.s);
			}
		} else {
			switch (l) {
			case -1:				//End of file
				break;
			case -2:
				//Length of read not equal to length of quality values
				nameLength = std::min(MAX_READNAME_LENGTH - 1, kseq->name.l);
				memcpy(read->name, kseq->name.s, nameLength);
				read->name[nameLength] = '\0';
				break;
			default:
				//Unknown error. Should not happen.
				throw "Unknown error while parsing. Please check whether the input file is corrupted!";
			}
		}
		return l;
	}

//	int copyToRead(SAMRecord * read, kseq_t * kseq, int const l) {
//
//		if (l >= 0) {
//			if (kseq->seq.l == kseq->qual.l || kseq->qual.l == 0) {
//				read->set_read_name(string(kseq->name.s));
//				read->set_sequence(string(kseq->seq.s));
//				read->set_qualities(string(kseq->qual.s));
//			} else {
//				throw "Error while parsing. Read length not equal to length of quality values!";
//			}
//		} else {
//			switch (l) {
//			case -1:				//End of file
//				break;
//			case -2:
//				//Length of read not equal to length of quality values
//				read->set_read_name(string(kseq->name.s));
//				break;
//			default:
//				//Unknown error. Should not happen.
//				throw "Unknown error while parsing. Please check whether the input file is corrupted!";
//			}
//		}
//		return l;
//	}

}
;

#endif /* PARSER_H_ */
