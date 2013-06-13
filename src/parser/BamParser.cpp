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
			memcpy(read->seq.s, al->QueryBases.c_str(), read->seq.l * sizeof(char));

			if (!al->Qualities.empty()) {
				//copy the qualities
				read->qual.l = al->Qualities.size();
				memcpy(read->qual.s, al->Qualities.c_str(), read->qual.l * sizeof(char));
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
