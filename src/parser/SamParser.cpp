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

void SamParser::init(char const * fileName, bool const keepTags) {
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
	tmp = kseq_init(fp);
	parseAdditionalInfo = keepTags;
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
int SamParser::doParseRead(MappedRead * read) {
	tmp->name.l = 0;
	tmp->seq.l = 0;
	tmp->qual.l = 0;

	while (gzgets(fp, buffer, buffer_size) != NULL) {
		char * lineBuffer = buffer;

		if (*lineBuffer != '@' && *lineBuffer != '\n') {
			//Name
			lineBuffer = readField(lineBuffer, tmp->name);
			if (*lineBuffer == '\0')
				return 0;

			//Skip one \t
			lineBuffer += 1;

			//Flags
			bool reverse = atoi(lineBuffer) & 0x10;

			int skip = 8;
			while (skip > 0 && *lineBuffer != '\0') {
				if (*lineBuffer++ == '\t')
					skip -= 1;
			}
			if (*lineBuffer == '\0')
				return 0;

			//Sequence
			lineBuffer = readField(lineBuffer, tmp->seq);
			if (reverse) {
				computeReverseSeq(tmp->seq.s, tmp->seq.l);
			}

			if (*lineBuffer == '\0')
				return 0;

			//Skip one \t
			lineBuffer += 1;

			//Quality
			lineBuffer = readField(lineBuffer, tmp->qual);
			if (reverse) {
				std::reverse(tmp->qual.s, &tmp->qual.s[strlen(tmp->qual.s)]);
			}

			if (tmp->qual.l == tmp->seq.l
					|| (tmp->qual.l == 1 && tmp->qual.s[0] == '*')) {

				if (parseAdditionalInfo) {
					int addInfoLen = strlen(lineBuffer) - 1;
					if (addInfoLen > 2) {
						if (read->AdditionalInfo == 0) {
							read->AdditionalInfo = new char[addInfoLen + 1];
							memcpy(read->AdditionalInfo, lineBuffer,
									addInfoLen);
							read->AdditionalInfo[addInfoLen] = '\0';
						}
					}
				}
				return copyToRead(read, tmp, tmp->seq.l);
			} else {
				//if (tmp->qual.l == 1 && tmp->qual.s[0] == '*') {
				//	tmp->qual.l = 0;
				//	return 0;
				//} else {
				return copyToRead(read, tmp, -2);
				//}
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

int SamParser::doParseRead(SAMRecord * read) {

	while (gzgets(fp, buffer, buffer_size) != NULL) {
		char * lineBuffer = buffer;

		if (*lineBuffer != '@' && *lineBuffer != '\n') {
			//Name
			string name;
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				name += *lineBuffer;
				lineBuffer++;
			}
			read->set_read_name(name);
			name.clear();

			//Skip one \t
			lineBuffer ++;

			//Flags
			read->set_mapped_flag(atoi(lineBuffer));
		}

		if(read->is_mapped()){
			//Skip the rest:
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				lineBuffer++;
			}
			//Skip one \t
			lineBuffer ++;
			string chr;
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				chr += *lineBuffer;
				lineBuffer++;
			}
			read->set_chr(chr);
			chr.clear();
			//Skip one \t
			lineBuffer ++;
			read->set_mapping_pos(atoi(lineBuffer));
			//Skip the rest:
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				lineBuffer++;
			}
			//Skip one \t
			lineBuffer++;
			read->set_mapping_quality(atoi(lineBuffer));
			//Skip the rest:
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				lineBuffer++;
			}
			//Skip one \t
			lineBuffer ++;
			string cigar;
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				cigar += *lineBuffer;
				lineBuffer++;
			}
			read->set_CIGAR(cigar);
			cigar.clear();

			int skip = 4;
			while (skip > 0 && *lineBuffer != '\0') {
				if (*lineBuffer++ == '\t')
					skip -= 1;
			}
			if (*lineBuffer == '\0') {
				return 0;
			}
			string seq;
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				seq += *lineBuffer;
				lineBuffer++;
			}
			read->set_sequence(seq);
			seq.clear();
			//Skip one \t
			lineBuffer ++;
			string qual;
			while (*lineBuffer != '\0' && *lineBuffer != '\t') {
				qual += *lineBuffer;
				lineBuffer++;
			}
			read->set_qualities(qual);
			qual.clear();
			//Skip one \t
			lineBuffer++;
			string tags;
			while (*lineBuffer != '\0') {
				tags += *lineBuffer;
				lineBuffer++;
			}
			read->set_tags(tags);
			tags.clear();
			return read->get_sequence().size();
		}
	}
	return -1;
}

