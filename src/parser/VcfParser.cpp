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
	int len = 0;
	for (int i = 0; i < SequenceProvider.GetRefCount();i++) {
		refmap[std::string(SequenceProvider.GetRefName(i,len),len)]=SequenceProvider.GetRefStart(i);
	}
}

VcfParser::~VcfParser() {

}

uint VcfParser::getRefStart(std::string ref) {
	for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
		int len = 0;
		std::string test = SequenceProvider.GetRefName(i * 2, len);
		if ( test == ref) {
			return SequenceProvider.GetRefStart(i * 2);
		}
	}
	return 0;
}

void VcfParser::open(char const * fileName)
{
	static const uint buffer_size = 512;
	char buffer[ buffer_size ];
	memset(buffer, 0, sizeof(char) * buffer_size);

	gzFile fp = gzopen(fileName, "r");
	if (!fp) {
		Log.Error("Failed to open VCF file ",fileName);
		return;
	}

	std::string vcf_data = "";

	while (gzgets(fp, buffer, buffer_size) != NULL) {
		vcf_data += std::string(buffer);
		memset(buffer, 0, sizeof(char) * buffer_size);
	}

	std::string line_buffer = "";
	uint line_num = 1;

	for(uint i = 0; i < vcf_data.size(); ++ i ) {
		if(vcf_data[i]=='\n')
		{
			parse_line(line_buffer, line_num);
			line_buffer = "";
			line_num++;
		} else{
			line_buffer.push_back(vcf_data[i]);
		}
	}
	
	if( line_buffer != "" )
		parse_line(line_buffer, line_num);
}

void VcfParser::parse_line(std::string line, uint line_num)
{
	if(line.size()==0) return;
	if(line[0]=='#')   return;

	//Trim
	uint begin = 0;
	uint end = line.size();

	for(uint i = 0; i < line.size(); ++ i )
	{
		if( line[i] == '\t' || line[i] == '\r' || line[i] == ' ' )
			begin ++;
		else
			break;
	}

	for(int i = end - 1; i >= 0; -- i )
	{
		if( line[i] == '\t' || line[i] == '\r' || line[i] == ' ' )
			end --;
		else
			break;
	}

	line = line.substr(begin,end-begin);

	//Split
	std::string buffer;
	std::vector<std::string> parts;

	for(uint i = 0; i < line.size(); ++ i )
	{
		if( line[i] == '\t' )
		{
			if( buffer != "" )
				parts.push_back(buffer);
			else
				Log.Message("Unexpected tab delimiter in VCF file, line %u", line_num);

			buffer = "";
		} else {
			buffer.push_back(line[i]);
		}
	}

	if( buffer != "" )
		parts.push_back(buffer);

	if( parts.size() < 8 )
	{
		Log.Message("Field count < 8 in VCF file, line %u", line_num);
		return;
	}

	std::string chrom = parts[0];
	std::string pos = parts[1];
	std::string ref = parts[3];
	std::string alt = parts[4];
	std::string current_alt = "";

	if( ref == "." )
	{
		//TODO: Needed?
		//return;
	}

	for( uint i = 0; i < parts[4].size(); ++ i )
	{
		if(alt[i]== ',')
		{
			if( current_alt != "" )
				add_line(chrom,pos,ref,current_alt,line_num);
			current_alt = "";
		} else {
			current_alt.push_back(alt[i]);
		}
	}

	if( current_alt != "" )
		add_line(chrom,pos,ref,current_alt,line_num);
}

void VcfParser::add_line(std::string chrom, std::string pos, std::string ref, std::string alt, uint line_num)
{
	VcfSNP snp;

	if(refmap.find(chrom) == refmap.end() )
	{
		Log.Message("Chromosome '%s' not found in reference but in VCF file, line %u", chrom.c_str(), line_num );
		return;
	}

	if( alt == "." )
	{
		//Missing alternative
		return;
	}

	if( ! isSequence( ref ) || ! isSequence( alt ) )
	{
		return;
	}

	snp.pos = getRefStart(chrom) + atoi(pos.c_str());
	snp.alt = alt;
	snp.ref = ref;
	snps.push_back(snp);
}

bool VcfParser::isSequence(const std::string& what)
{
	for(uint i = 0; i < what.size(); ++ i )
	{
		switch(what[i])
		{
			case 'A':
			case 'C':
			case 'T':
			case 'G':
			case 'N':
			break;

			default: return false;
		}
	}

	return true;
}