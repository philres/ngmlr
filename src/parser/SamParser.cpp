/*
 * SamParser.cpp
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#include "SamParser.h"

#include <algorithm>

#include <string.h>

#include "Log.h"

void SamParser::init(char const * fileName) {
	fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("File does not exist ",fileName);
	}

	buffer = new char[buffer_size];

	//parse_all = Config.Exists("parse_all") && Config.GetInt("parse_all") == 1;
	parse_all = true;
	if (!parse_all) {
		Log.Warning("Skipping all mapped reads in SAM file.");
	}
	read = kseq_init(fp);
}

static inline char * readField(char * lineBuffer, kstring_t & str) {
	while (*lineBuffer != '\t' && *lineBuffer != '\0' && *lineBuffer != '\n') {
		if (str.l + 1 > str.m) {
			str.m = str.l + 2;
			kroundup32(str.m);
			str.s = (char*) realloc(str.s, str.m);
		}
		str.s[str.l++] = *lineBuffer++;
	}
	str.s[str.l] = '\0';
	return lineBuffer;
}

static inline char cpl(char c) {
	if (c == 'A')
		return 'T';
	else if (c == 'T')
		return 'A';
	else if (c == 'C')
		return 'G';
	else if (c == 'G')
		return 'C';
	else
		return c;
}

//// swaps two bases and complements them
//static inline void rc(char & c1, char & c2)
//{
//	char x = c1;
//	c1 = cpl(c2);
//	c2 = cpl(x);
//}

void computeReverseSeq(char * Seq, int qryMaxLen) {
	char * RevSeq = new char[qryMaxLen + 1];
	memset(RevSeq, 0, qryMaxLen + 1);
	memcpy(RevSeq, Seq, qryMaxLen);

	char * fwd = RevSeq;
	char * rev = Seq + qryMaxLen - 1;

	for (int i = 0; i < qryMaxLen; ++i) {
		*rev-- = cpl(*fwd++);
	}
	//if()
	delete RevSeq;
	RevSeq = 0;
}

/* Return value:
 >=0  length of the sequence (normal)
 -1   end-of-file
 -2   truncated quality string
 */
size_t SamParser::parseRead() {
	read->name.l = 0;
	read->seq.l = 0;
	read->qual.l = 0;

	while (gzgets(fp, buffer, buffer_size) != NULL) {
		char * lineBuffer = buffer;

		//Log.Message("%s", lineBuffer);
		if (*lineBuffer != '@' && *lineBuffer != '\n') {

			//Name
			lineBuffer = readField(lineBuffer, read->name);
			if (*lineBuffer == '\0')
				return 0;

			//Skip one \t
			lineBuffer += 1;

			//Flags
			bool reverse = atoi(lineBuffer) & 0x10;

			//Log.Message("%d %d", atoi(lineBuffer), reverse);
			int skip = 8;
			while (skip > 0 && *lineBuffer != '\0') {
				if (*lineBuffer++ == '\t')
					skip -= 1;
			}
			if (*lineBuffer == '\0')
				return 0;

			//Sequence
			lineBuffer = readField(lineBuffer, read->seq);
			if (reverse) {
				//char * tmp = read->seq.s;
				//Log.Message("Seq:    %s", read->seq.s);
				computeReverseSeq(read->seq.s, read->seq.l);
				//Log.Message("RevSeq: %s", read->seq.s);

				//delete tmp; tmp = 0;
			}

			if (*lineBuffer == '\0')
				return 0;

			//Skip one \t
			lineBuffer += 1;

			//Quality
			lineBuffer = readField(lineBuffer, read->qual);
			if (reverse) {
				std::reverse(read->qual.s, &read->qual.s[strlen(read->qual.s)]);
			}

			if (read->qual.l == read->seq.l) {
				return read->seq.l;
			} else {
				if (read->qual.l == 1 && read->qual.s[0] == '*') {
					read->qual.l = 0;
					return 0;
				} else {
					return -2;
				}
			}
		}
	}
	return -1;
}

///* Return value:
// >=0  length of the sequence (normal)
// -1   end-of-file
// -2   truncated quality string
// */
//size_t SamParser::parseRead(kseq_t *& read) {
//	while (gzgets(fp, buffer, buffer_size) != NULL) {
//		read->name.s = &buffer[0]; //name
//		int count = 0;
//		int len = 0;
//		for (int i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
//			if (count == 1 && buffer[i - 1] == '\t') {
//				std::string a;
//				a = buffer[i];
//				read->name.l = i - 2;
//				read->name.m = i - 2;
//				if (!(atoi(&buffer[i]) & 0x4) && !parse_all) {
//					break; //TODO print line to output file
//				}
//			}
//			if (count == 9 && buffer[i - 1] == '\t') {
//
//				read->seq.s = &buffer[i];
//			}
//			if (count == 9 && buffer[i] != '\t') {
//				len++;
//			}
//			if (count == 10 && buffer[i - 1] == '\t') { //point to the qv values;
//				read->seq.l = len;
//				read->seq.m = len;
//				if (buffer[i] != '*') {
//					read->qual.s = &buffer[i];
//					read->qual.l = len;
//					read->qual.m = len;
//				} else {
//					read->qual.l = 0;
//					read->qual.m = 0;
//				}
//				return len;
//			}
//			if (buffer[i] == '\t') {
//				count++;
//			}
//		}
//	}
//	return -1;
//}

