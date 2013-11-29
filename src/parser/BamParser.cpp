/*
 * BamParser.cpp
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#include "BamParser.h"

#ifdef _BAM
#include "Config.h"
#include "Log.h"
#include <stdio.h>
#include <string.h>

void BamParser::init(char const * fileName) {
	std::vector<std::string> tmps;
	tmps.push_back(fileName);

	if (!reader.Open(tmps)) {
		Log.Error("File does not exist ",fileName);
	}
	al = new BamAlignment();
	parse_all = bool(Config.Exists("parse_all") && Config.GetInt("parse_all") == 1);

	read = kseq_init(fp);
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

/* Return value:
 >=0  length of the sequence (normal)
 -1   end-of-file
 -2   truncated quality string
 */
size_t BamParser::parseRead() {

	while (reader.GetNextAlignmentCore(al[0])) {
		if (!al->IsMapped() || parse_all) {
			al->BuildCharData();
			if (al->Name.size() > read->name.l) { //parse the name
				read->name.m = al->Name.size();
				kroundup32(read->name.m);
				// round to the next k^2
				read->name.s = (char*) realloc(read->name.s, read->name.m);
			}
			//copy the name
			read->name.l = al->Name.size();
			memcpy(read->name.s, al->Name.c_str(), read->name.l * sizeof(char));
			if (al->QueryBases.size() > read->seq.m) { //adjust the sequence size
				//m is size of read
				read->seq.m = al->QueryBases.size();
				kroundup32(read->seq.m);
				// round to the next k^2
				read->seq.s = (char*) realloc(read->seq.s, read->seq.m);

				if (!al->Qualities.empty()) {
					read->qual.m = al->Qualities.size();
					kroundup32(read->qual.m);
					read->qual.s = (char*) realloc(read->qual.s, read->qual.m);
				}
			}
			//copy the sequence
			read->seq.l = al->QueryBases.size();
			if(al->IsReverseStrand()) {
				char const * fwd = al->QueryBases.c_str();
				char * rev = read->seq.s + read->seq.l - 1;

				for (int i = 0; i < read->seq.l; ++i) {
					*rev-- = cpl(*fwd++);
				}
			} else {
				memcpy(read->seq.s, al->QueryBases.c_str(), read->seq.l * sizeof(char));
			}

			if (!al->Qualities.empty()) {
				//copy the qualities
				read->qual.l = al->Qualities.size();
				if(al->IsReverseStrand()) {
					for(int i = 0; i < read->qual.l; ++i) {
						read->qual.s[i] = al->Qualities.c_str()[read->qual.l - 1 - i];
					}
				} else {
					memcpy(read->qual.s, al->Qualities.c_str(), read->qual.l * sizeof(char));
				}
			}

			if (read->seq.l != read->qual.l) {
				return -2;
			}
			return read->seq.m;
		} else {

			//TODO: print directly into the output files

		}
	}
	return -1;
}

#endif
