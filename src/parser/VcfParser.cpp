/*
 * VcfParser.h
 *
 * Created originally on: Sep 11, 2012 by Fritz
 * Rewrite: December 19, 2014
 *      Author: moritz
 */

#include "VcfParser.h"
#include "SequenceProvider.h"
#include "Config.h"
#include "Log.h"
#include <zlib.h>

VcfParser::VcfParser() {

}

VcfParser::~VcfParser() {

}

void VcfParser::open(char const * fileName)
{
	static const uint buffer_size = 512;
	char buffer[ buffer_size ];

	gzFile fp = gzopen(fileName, "r");
	if (!fp) {
		//File does not exist
		Log.Error("Failed to open VCF file ",fileName);
		return;
	}

	
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
				snps.push_back(snp);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
	}
	refmap.clear();	
}