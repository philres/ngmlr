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

void BamParser::init(char const * fileName, bool const keepTags) {
	std::vector<std::string> tmps;
	tmps.push_back(fileName);

	if (!reader.Open(tmps)) {
		Log.Error("File does not exist ",fileName);
	}
	al = new BamAlignment();
	parse_all = bool(
			Config.Exists("parse_all") && Config.GetInt("parse_all") == 1);

	tmp = kseq_init(fp);
	parseAdditionalInfo = keepTags;
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
int BamParser::doParseRead(MappedRead * read) {

	while (reader.GetNextAlignmentCore(al[0])) {
		if (!al->IsMapped() || parse_all) {
			al->BuildCharData();
			if (al->Name.size() > tmp->name.l) { //parse the name
				tmp->name.m = al->Name.size();
				kroundup32(tmp->name.m);
				// round to the next k^2
				tmp->name.s = (char*) realloc(tmp->name.s, tmp->name.m);
			}
			//copy the name
			tmp->name.l = al->Name.size();
			memcpy(tmp->name.s, al->Name.c_str(), tmp->name.l * sizeof(char));
			if (al->QueryBases.size() > tmp->seq.m) { //adjust the sequence size
				//m is size of read
				tmp->seq.m = al->QueryBases.size();
				kroundup32(tmp->seq.m);
				// round to the next k^2
				tmp->seq.s = (char*) realloc(tmp->seq.s, tmp->seq.m);

				if (!al->Qualities.empty()) {
					tmp->qual.m = al->Qualities.size();
					kroundup32(tmp->qual.m);
					tmp->qual.s = (char*) realloc(tmp->qual.s, tmp->qual.m);
				}
			}
			//copy the sequence
			tmp->seq.l = al->QueryBases.size();
			if (al->IsReverseStrand()) {
				char const * fwd = al->QueryBases.c_str();
				char * rev = tmp->seq.s + tmp->seq.l - 1;

				for (size_t i = 0; i < tmp->seq.l; ++i) {
					*rev-- = cpl(*fwd++);
				}
			} else {
				memcpy(tmp->seq.s, al->QueryBases.c_str(),
						tmp->seq.l * sizeof(char));
			}

			if (!al->Qualities.empty()) {
				//copy the qualities
				tmp->qual.l = al->Qualities.size();
				if (al->IsReverseStrand()) {
					for (size_t i = 0; i < tmp->qual.l; ++i) {
						tmp->qual.s[i] = al->Qualities.c_str()[tmp->qual.l - 1
								- i];
					}
				} else {
					memcpy(tmp->qual.s, al->Qualities.c_str(),
							tmp->qual.l * sizeof(char));
				}
			}

			if (tmp->qual.l == tmp->seq.l
					|| (tmp->qual.l == 1 && tmp->qual.s[0] == '*')) {

				if (parseAdditionalInfo) {
					size_t position = 0;

					std::vector<std::string> tags = al->GetTagNames();

					for (size_t i = 0; i < tags.size(); i++) {

						char type = 0;
						al->GetTagType(tags[i], type);
						if (type == Constants::BAM_TAG_TYPE_INT8
								|| type == Constants::BAM_TAG_TYPE_INT16
								|| type == Constants::BAM_TAG_TYPE_INT32) {
							int value = 0;
							al->GetTag<int>(tags[i], value);
							position += sprintf(additionalInfo + position,
									"\t%s:%c:%d", tags[i].c_str(), type, value);
						} else if (type == Constants::BAM_TAG_TYPE_UINT8
								|| type == Constants::BAM_TAG_TYPE_UINT16
								|| type == Constants::BAM_TAG_TYPE_UINT32) {
							uint value = 0;
							al->GetTag<uint>(tags[i], value);
							position += sprintf(additionalInfo + position,
									"\t%s:%c:%u", tags[i].c_str(), type, value);
						} else if (type == Constants::BAM_TAG_TYPE_STRING
								|| type == Constants::BAM_TAG_TYPE_ASCII) {
							std::string value;
							al->GetTag<std::string>(tags[i], value);
							position += sprintf(additionalInfo + position,
									"\t%s:%c:%s", tags[i].c_str(), type,
									value.c_str());
						} else if (type == Constants::BAM_TAG_TYPE_FLOAT) {
							float value = 0;
							al->GetTag<float>(tags[i], value);
							position += sprintf(additionalInfo + position,
									"\t%s:%c:%f", tags[i].c_str(), type, value);
						}

//					switch ( type ) {
//						case (Constants::BAM_TAG_TYPE_HEX) :
//						break;
//						case (Constants::BAM_TAG_TYPE_ARRAY) :
//						break;
//					}

					}
					if (position > 0) {
						read->AdditionalInfo = new char[position + 1];
						memcpy(read->AdditionalInfo, additionalInfo, position);
						read->AdditionalInfo[position] = '\0';
					}
				}
				return copyToRead(read, tmp, tmp->seq.l);
			} else {
				return copyToRead(read, tmp, -2);
			}

		} else {

			//TODO: print directly into the output files

		}
	}
	return -1;
}



int BamParser::doParseRead(SAMRecord * read) {

	while (reader.GetNextAlignmentCore(al[0])) {
		if (al->IsMapped()) {
			al->BuildCharData();
			read->set_read_name(al->Name);
			read->set_sequence( al->QueryBases);
			read->set_CIGAR(al->CigarData);
			read->set_chr(reader.GetReferenceData()[al->RefID].RefName);
			read->set_mapped_flag(al->AlignmentFlag);
			read->set_mapping_pos(al->Position);
			read->set_mapping_quality(al->MapQuality);
			read->set_qualities(al->Qualities);
			read->set_sequence(al->QueryBases);
			read->set_tags(al->TagData);
			return al->QueryBases.size();
		}
	}

	return -1;
}

#endif
