/*
 * filter.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: philipp_
 */

#include "filter.h"

#include <cstring>
#include <cmath>

#include "SamFile.h"
#include "SamValidation.h"

#include "Log.h"

using std::string;

//first: read (raw)
//second: ref (wird erzeugt)

//void Alignment::computeAlignment() {
//	int pos = 0;
//	bool flag = false;
//
//	for (size_t i = 0; i < al->CigarData.size(); i++) {
//		if (al->CigarData[i].Type == 'I') {
//			for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
//				alignment.second.insert(pos, "-");
//				alignment.second.erase(alignment.second.size() - 1, 1);
//				pos++;
//			}
//		} else if (al->CigarData[i].Type == 'D') {
//			for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
//				alignment.first.insert(pos, "-");
//				pos++;
//			}
//		} else if (al->CigarData[i].Type == 'S') {
//
//			if (pos == 0) { //front side
//				alignment.second.erase(
//						((int) alignment.second.size())
//								- al->CigarData[i].Length,
//						al->CigarData[i].Length);
//			} else { //backside
//				alignment.second.erase(pos, al->CigarData[i].Length);
//			}
//			alignment.first.erase(pos, al->CigarData[i].Length);
//
//		} else if (al->CigarData[i].Type == 'M') {
//
//			flag = true;
//			pos += al->CigarData[i].Length;
//		} else if (al->CigarData[i].Type == 'H') {
//
//		}
//	}
//	for (size_t i = 0; i < alignment.first.size(); i++) {
//		if (alignment.first[i] == '=') {
//			alignment.first[i] = alignment.second[i];
//		}
//	}
//	is_computed = true;
//
//	if (alignment.first.size() != alignment.second.size()) {
//		cerr << "Error alignment has different length" << endl;
//		cerr << " ignoring alignment " << al->Name << endl;
//		cerr << al->Position << endl;
//
//		cerr << endl;
//		cerr << "read: " << alignment.first << endl;
//		cerr << endl;
//		cerr << " ref: " << alignment.second << endl;
//		cerr << endl;
//		cerr << orig_length << endl;
//		vector < CigarOp > cig = getCigar();
//
//		for (size_t i = 0; i < cig.size(); i++) {
//			cerr << cig[i].Length << cig[i].Type << " ";
//		}
//		cerr << endl;
//		exit(0);
//		return;
//	}
//
//}

int filter(int argc, char **argv) {

	SamFile samIn;
	samIn.OpenForRead(argv[1]);

	SamStatus::Status returnStatus = SamStatus::SUCCESS;

	// Read the sam header.
	SamFileHeader samHeader;
	samIn.ReadHeader(samHeader);

	SamFile samOut;
	samOut.OpenForWrite("output.sam");
	samOut.WriteHeader(samHeader);

	SamRecord samRecord;
	while (samIn.ReadRecord(samHeader, samRecord)) {
		Log.Message("Name: %s", samRecord.getReadName());
		Log.Message("Seq:  %s", samRecord.getSequence());
		Log.Message("NM:   %d", samRecord.getInteger("NM"));
		Log.Message("Cigar:%s", samRecord.getCigar());
		Log.Message("MD:   %s", samRecord.getString("MD").c_str());

		Cigar * cigar = samRecord.getCigarInfo();
		string cigarString;
		cigar->getExpandedString(cigarString);

		int alignmentLength = cigarString.length();

		//Expand MD
		char * md = new char[alignmentLength + 1];
		memset(md, '=', alignmentLength);
		md[alignmentLength] = '\0';

		const char * mdString = samRecord.getString("MD").c_str();
		int mdLength = samRecord.getString("MD").Length();



		Log.Message("%d",  samRecord.get0BasedPosition() - samRecord.get0BasedUnclippedStart());
		int i = 0;
		int mdIndex = samRecord.get0BasedPosition() - samRecord.get0BasedUnclippedStart();
		while(i < mdLength) {
			if(mdString[i] == '^') {
				i += 1;
			} else if(mdString[i] >= '0' && mdString[i] <= '9') {
				int offset = atoi(mdString + i);
				for(int j = 0; j < offset; ++j) {
					//Log.Green("cigar at %d is %c", mdIndex, cigarString[mdIndex]);
					if(cigarString[mdIndex] == 'I') {
						Log.Green("jetzt");
						md[mdIndex++] = 'I';
					}
					mdIndex += 1;
				}
				//mdIndex += offset;
				if(offset > 0) {
					i += (int)log10(offset) + 1;
				} else {
					i += 1;
				}

			} else {
				while(cigarString[mdIndex] == 'I') {
					md[mdIndex++] = 'I';
				}
				//Log.Message("Setting %d to %c", mdIndex, mdString[i]);
				md[mdIndex++] = mdString[i++];
			}
		}



		Log.Message("Cigar:%s", cigarString.c_str());
		Log.Message("MD:   %s", md);

		char * ref = new char[alignmentLength + 1];
		char * read = new char[alignmentLength + 1];
		memset(ref, '-', alignmentLength);
		memset(read, '-', alignmentLength);
		ref[alignmentLength] = '\0';
		read[alignmentLength] = '\0';

		int refIndex = 0;
		int readIndex = 0;

		for(int i = 0; i < cigarString.length(); ++i) {
			switch(cigarString[i]) {
				case 'M':
				ref[i] = read[i] = samRecord.getSequence(readIndex++);
				if(md[i] != '=' && md[i] != 'I') {
					ref[i] = md[i];
				}
				break;
				case 'D':
				ref[i] = md[i];
				break;
				case 'I':
				read[i] = samRecord.getSequence(readIndex++);
				break;
			}
		}
		ref[cigarString.length()] = '\0';
		read[cigarString.length()] = '\0';

		delete[] md; md = 0;

		Log.Message("Ref:  %s", ref);
		Log.Message("Read: %s", read);

		delete[] ref; ref = 0;
		delete[] read; read = 0;
		getchar();
//if(samRecord.getInteger("NM") <= 2) {
//	samOut.WriteRecord(samHeader, samRecord);
//}
	}

}
