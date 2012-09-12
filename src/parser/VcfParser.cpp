/*
 * VcfParser.cpp
 *
 *  Created on: Sep 11, 2012
 *      Author: fritz
 */

#include "VcfParser.h"
#include "Config.h"
#include "Log.h"
#include "SequenceProvider.h"
#include <list>



bool VcfParser::compare_SNP(SNP first, SNP second){
	return first.pos< second.pos;
}
void VcfParser::VcfParser(char const * fileName) {
	fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("File does not exist ",fileName);
	}

	buffer = new char[buffer_size];
	std::map<std::string, int> refmap;

	int len = 0;
	for (int i = 0; i < SequenceProvider.GetRefCount();i++) {
		refmap[std::string(SequenceProvider.GetRefName(i,len),len)]=SequenceProvider.GetRefStart(i);
	}

	//run over header:
	while (gzgets(fp, buffer, buffer_size) > 0 && buffer[0] == '#') {
	}



	while (gzgets(fp, buffer, buffer_size) != NULL) {

		int count = 0;
		std::string name = "";
		SNP snp;
		for (int i = 0;	i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';i++) {

			if (count == 0 && buffer[i] != '\t') {
				name += buffer[i];
			}

			if (count == 2 && buffer[i - 1] == '\t') {
				snp.pos = atoi(&buffer[i]);
				snp.pos += refmap[name.c_str()]; //adjust for concatenated Refs
			}
			if (count == 4 && buffer[i] != '\t') {
				snp.alt = buffer[i];
				snp_list.push_back(snp);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
	}
	refmap.clear();
	snp_list.sort(compare_SNP);

}

SNP VcfParser::getNextSnp() {
	SNP tmp = snp_list.front();
	snp_list.pop_front();
	return tmp;
}
