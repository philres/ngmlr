/*
 * SamParser.cpp
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#include "SamParser.h"
#include "Config.h"
#include "Log.h"

kseq_t * SamParser::init_seq(char const * fileName) {
	fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("File does not exist ",fileName);
	}

	buffer = new char[buffer_size];
	//run over header:
	while (gzgets(fp, buffer, buffer_size) > 0 && buffer[0] == '@') {
	}

	parse_all = Config.Exists("parse_all") && Config.GetInt("parse_all") == 1;
	if (!parse_all) {
		Log.Warning("Skipping all mapped reads in SAM file.");
	}
	kseq_t *seq = kseq_init(fp);
	return seq;
}

/* Return value:
 >=0  length of the sequence (normal)
 -1   end-of-file
 -2   truncated quality string
 */
size_t SamParser::parseRead(kseq_t *& read) {
	while (gzgets(fp, buffer, buffer_size) != NULL) {
		read->name.s = &buffer[0]; //name
		int count = 0;
		int len = 0;
		for (int i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 1 && buffer[i - 1] == '\t') {
				std::string a;
				a = buffer[i];
				read->name.l = i - 2;
				read->name.m = i - 2;
				if (!(atoi(&buffer[i]) & 0x4) && !parse_all) {
					break; //TODO print line to output file
				}
			}
			if (count == 9 && buffer[i - 1] == '\t') {
				read->seq.s = &buffer[i];
			}
			if (count == 9 && buffer[i] != '\t') {
				len++;
			}
			if (count == 10 && buffer[i - 1] == '\t') { //point to the qv values;
				read->seq.l = len;
				read->seq.m = len;
				if (buffer[i] != '*') {
					read->qual.s = &buffer[i];
					read->qual.l = len;
					read->qual.m = len;
				} else {
					read->qual.l = 0;
					read->qual.m = 0;
				}
				return len;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
	}
	return -1;
}

