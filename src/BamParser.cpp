/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
 */

#include "BamParser.h"

#include "IConfig.h"
#include "Log.h"
#include <stdio.h>
#include <string.h>

#include <string>
#include <iostream>
#include <fstream>

#undef module_name
#define module_name "BAM"

bool isPrimary(BamAlignment * al) {
	return !(al->AlignmentFlag & 0x100) && !(al->AlignmentFlag & 0x800);
}

void BamParser::parseBedFile(char const * fileName) {
	Log.Message("Parsing BED file: %s", fileName);

	std::ifstream input(fileName);
	std::string line;

	while (std::getline(input, line)) {
		if(line.length() > 5) {

			BamRegion region;

			std::string delimiter = "\t";
			int pos = line.find(delimiter);
			std::string token = line.substr(0, pos); // token is "scott"

			if(reader.GetReferenceID(token) < 0) {
				Log.Error("Invalid chromosome name!");
			}
			region.LeftRefID = reader.GetReferenceID(token);
			region.RightRefID = reader.GetReferenceID(token);

			int lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos); // token is "scott"

			region.LeftPosition = atoi(token.c_str());

			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			region.RightPosition = atoi(token.c_str());


			regions.push_back(region);

			while(pos != std::string::npos) {
				lastPos = pos + 1;
				pos = line.find(delimiter, lastPos);

				token = line.substr(lastPos, pos - lastPos); // token is "scott"
				readNames[token] = 1;
				Log.Message("Adding name: %s", token.c_str());
			}

		}
	}
}

void BamParser::parseRealignFile(char const * fileName) {
	Log.Message("Parsing Realign file: %s", fileName);

	std::ifstream input(fileName);
	std::string line;

	while (std::getline(input, line)) {
		if(line.length() > 5) {

			int pos = 0;
			int lastPos = 0;

			BamRegion region;

			std::string delimiter = "\t";

			pos = line.find(delimiter);
			std::string token = line.substr(0, pos); // token is "scott"


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			if(reader.GetReferenceID(token) < 0) {
				Log.Error("Invalid chromosome name!");
			}
			region.LeftRefID = reader.GetReferenceID(token);
			region.RightRefID = region.LeftRefID;


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos); // token is "scott"

			region.LeftPosition = atoi(token.c_str());
			region.RightPosition = region.LeftPosition + 1;

			regions.push_back(region);

			BamRegion region2;


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);

			region2.LeftRefID = reader.GetReferenceID(token);
			region2.RightRefID = region.LeftRefID;


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			region2.LeftPosition = atoi(token.c_str());
			region2.RightPosition = region.LeftPosition + 1;


			regions.push_back(region2);

			while(pos != std::string::npos) {
				lastPos = pos + 1;
				pos = line.find(delimiter, lastPos);

				token = line.substr(lastPos, pos - lastPos); // token is "scott"
				readNames[token] = 1;
				Log.Message("Adding name: %s", token.c_str());
			}

		}
	}
}

void BamParser::parseSnifflesFile(char const * fileName) {

	Log.Message("Parsing sniffles file: %s", fileName);

	std::ifstream input(fileName);
	std::string line;

	while (std::getline(input, line)) {
		if(line.length() > 5) {

			BamRegion region;

			std::string delimiter = "\t";
			int pos = line.find(delimiter);
			std::string token = line.substr(0, pos); // token is "scott"

			if(reader.GetReferenceID(token) < 0) {
				Log.Error("Invalid chromosome name!");
			}
			region.LeftRefID = reader.GetReferenceID(token);
			region.RightRefID = reader.GetReferenceID(token);


			int lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos); // token is "scott"

			region.LeftPosition = atoi(token.c_str());


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			region.RightPosition = atoi(token.c_str());


			regions.push_back(region);

			BamRegion region2;

			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			region2.LeftRefID = reader.GetReferenceID(token);
			region2.RightRefID = reader.GetReferenceID(token);


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			region2.LeftPosition = atoi(token.c_str());


			lastPos = pos + 1;
			pos = line.find(delimiter, lastPos);
			token = line.substr(lastPos, pos - lastPos);// token is "scott"

			region.RightPosition = atoi(token.c_str());


			regions.push_back(region2);

			while(pos != std::string::npos) {
				lastPos = pos + 1;
				pos = line.find(delimiter, lastPos);

				token = line.substr(lastPos, pos - lastPos); // token is "scott"
				readNames[token] = 1;
				Log.Message("Adding name: %s", token.c_str());
			}

		}
	}

}

void BamParser::init(char const * fileName, bool const keepTags) {
	std::vector<std::string> tmps;
	tmps.push_back(fileName);

	if (!reader.Open(tmps)) {
		Log.Error("File does not exist ",fileName);
	}
	al = new BamAlignment();

	reader.LocateIndexes();

	char const * const bedFile = Config.getBedFile();
	if (bedFile != 0) {
		parseBedFile(bedFile);
		if (regions.size() > 0) {
			BamRegion region = regions.back();

			Log.Message("Setting region to %d:%d-%d", region.LeftRefID, region.LeftPosition, region.RightPosition);
			reader.SetRegion(region);
			regions.pop_back();
		}
	}

	tmp = kseq_init(fp);
	parseAdditionalInfo = keepTags;

	Log.Message("BAM parser initialized");

//	int readOffset = Config.GetInt(READ_OFFSET);
//	if (readOffset > 0) {
//
//		Log.Message("Skipping first %d reads", readOffset);
//		int count = 0;
//		while ((reader.GetNextAlignmentCore(al[0])) && count < readOffset) {
//			if (isPrimary(al)) { // && readNames.find(al->Name) != readNames.end()) {
//				count += 1;
//			}
//		}
//		Log.Verbose("%d reads skipped", count);
//	}
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

int BamParser::doParseSingleRead(MappedRead * read, BamAlignment * al) {

	parsedReads += 1;
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
	tmp->name.s[tmp->name.l] = '\0';
	Log.Verbose("Read %s mapping to %d:%d parsed (flags: %d)", tmp->name.s, al->RefID, al->Position, al->AlignmentFlag);
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
		memcpy(tmp->seq.s, al->QueryBases.c_str(), tmp->seq.l * sizeof(char));
	}

	if (!al->Qualities.empty() && al->Qualities[0] != -1) {
		//copy the qualities
		tmp->qual.l = al->Qualities.size();
		if (al->IsReverseStrand()) {
			for (size_t i = 0; i < tmp->qual.l; ++i) {
				tmp->qual.s[i] = al->Qualities.c_str()[tmp->qual.l - 1 - i] + 0;
			}
		} else {
			memcpy(tmp->qual.s, al->Qualities.c_str(),
					tmp->qual.l * sizeof(char));
		}
	} else {
		tmp->qual.l = 0;
	}

	if (tmp->qual.l == tmp->seq.l
			|| (tmp->qual.l == 1 && tmp->qual.s[0] == '*') || tmp->qual.l == 0) {

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
					position += sprintf(additionalInfo + position, "\t%s:%c:%d",
							tags[i].c_str(), type, value);
				} else if (type == Constants::BAM_TAG_TYPE_UINT8
						|| type == Constants::BAM_TAG_TYPE_UINT16
						|| type == Constants::BAM_TAG_TYPE_UINT32) {
					uint value = 0;
					al->GetTag<uint>(tags[i], value);
					position += sprintf(additionalInfo + position, "\t%s:%c:%u",
							tags[i].c_str(), type, value);
				} else if (type == Constants::BAM_TAG_TYPE_STRING
						|| type == Constants::BAM_TAG_TYPE_ASCII) {
					std::string value;
					al->GetTag<std::string>(tags[i], value);
					position += sprintf(additionalInfo + position, "\t%s:%c:%s",
							tags[i].c_str(), type, value.c_str());
				} else if (type == Constants::BAM_TAG_TYPE_FLOAT) {
					float value = 0;
					al->GetTag<float>(tags[i], value);
					position += sprintf(additionalInfo + position, "\t%s:%c:%f",
							tags[i].c_str(), type, value);
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
}

/* Return value:
 >=0  length of the sequence (normal)
 -1   end-of-file
 -2   truncated quality string
 */
int BamParser::doParseRead(MappedRead * read) {

//	static int readNumber = Config.GetInt(READ_NUMBER);

	bool found = false;
	while (((found = reader.GetNextAlignmentCore(al[0])) || regions.size() > 0)) {
//			&& (readNumber == -1 || parsedReads < readNumber)) {
		if (!found) {
			BamRegion region = regions.back();
			Log.Message("Setting region to %d:%d-%d", region.LeftRefID, region.LeftPosition, region.RightPosition);
			reader.SetRegion(region);
			regions.pop_back();
		} else {
			if (isPrimary(al)) { // && readNames.find(al->Name) != readNames.end()) {
				al->BuildCharData();
				//					Log.Message("%d Looking for %s", x++, al->Name.c_str());
				if (readNames.size() == 0
						|| readNames.find(al->Name) != readNames.end()) {
					return doParseSingleRead(read, al);
				}
			}
		}
	}

	return -1;
}

int BamParser::doParseRead(SAMRecord * read) {

	while (reader.GetNextAlignmentCore(al[0])) {
		if (al->IsMapped()) {
			al->BuildCharData();
			read->set_read_name(al->Name);
			read->set_sequence(al->QueryBases);
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

