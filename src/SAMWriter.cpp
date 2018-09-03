/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "SAMWriter.h"

#include <string.h>
#include <sstream>
#include <algorithm>
#include <string.h>
#include <cmath>

#include "IConfig.h"
#include "SequenceProvider.h"
#include "NGM.h"
#include "Version.h"

//Format: http://samtools.sourceforge.net/SAM1.pdf
static int const report_offset = 1;

void SAMWriter::DoWriteProlog() {
	//TODO: check correct format

	char const * refName = 0;
	int refNameLength = 0;
	Print("@HD\tVN:1.0\tSO:unsorted\n");

	for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
		refName = SequenceProvider.GetRefName(i, refNameLength);
		Print("@SQ\tSN:%.*s\tLN:%llu\n", refNameLength, refName,SequenceProvider.GetRefLen(i));
		++i;
		m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer, false);
	}

	std::stringstream version;
	version << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_BUILD;
	Print("@PG\tID:ngmlr\tPN:nextgenmap-lr\tVN:%s\tCL:%s\n", version.str().c_str(), Config.getFullCommandLineCall());

	if (rgId != 0) {
		Print("@RG\tID:%s", rgId);

		char const * rgSm = Config.getRgSm();
		if(rgSm != 0)
			Print("\tSM:%s", rgSm);
		char const * rgLb = Config.getRgLb();
		if(rgLb != 0)
			Print("\tLB:%s", rgLb);
		char const * rgPl = Config.getRgPl();
		if(rgPl != 0)
			Print("\tPL:%s", rgPl);
		char const * rgDs = Config.getRgDs();
		if(rgDs != 0)
			Print("\tDS:%s", rgDs);
		char const * rgDt = Config.getRgDt();
		if(rgDt != 0)
			Print("\tDT:%s", rgDt);
		char const * rgPu = Config.getRgPu();
		if(rgPu != 0)
			Print("\tPU:%s", rgPu);
		char const * rgPi = Config.getRgPi();
		if(rgPi != 0)
			Print("\tPI:%s", rgPi);
		char const * rgPg = Config.getRgPg();
		if(rgPg != 0)
			Print("\tPG:%s", rgPg);
		char const * rgCn = Config.getRgCn();
		if(rgCn != 0)
			Print("\tCN:%s", rgCn);
		char const * rgFo = Config.getRgFo();
		if(rgFo != 0)
			Print("\tFO:%s", rgFo);
		char const * rgKs = Config.getRgKs();
		if(rgKs != 0)
			Print("\tKS:%s", rgKs);

		Print("\n");
	}

	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer, true);
}

void SAMWriter::DoWriteRead(MappedRead const * const read, int const scoreID) {
	DoWriteReadGeneric(read, scoreID, "*", -1, 0, read->mappingQlty);
	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer);
}

void SAMWriter::DoWriteReadGeneric(MappedRead const * const read,
		int const scoreID, char const * pRefName, int const pLoc,
		int const pDist, int const mappingQlty, int flags) {

	NGM.AddWrittenRead(read->ReadId);

	static bool const hardClip = Config.getHardClip();

	char const * readseq = read->Seq;
	char * readname = read->name;
	char * qltystr = read->qlty;

	//if (scoreID != 0) {
	if(!read->Alignments[scoreID].primary) {
		flags |= 0x800;
	}

	if (read->Scores[scoreID].Location.isReverse()) {
		readseq = read->RevSeq;
		if (qltystr != 0 && strlen(qltystr) > 0) {
			std::reverse(qltystr, &qltystr[read->length]);
		}
		flags |= 0x10;
	}

	int refnamelen = 0;
	char const * refname = SequenceProvider.GetRefName(read->Scores[scoreID].Location.getrefId(), refnamelen);

	//mandatory fields
	Print("%s\t", readname);
	Print("%d\t", flags);
	Print("%.*s\t", refnamelen, refname);
	Print("%u\t", read->Scores[scoreID].Location.m_Location + report_offset );
//	Print("%d\t", mappingQlty);
	Print("%d\t", read->Alignments[scoreID].MQ);

	int printLongCigar = (Config.getBamCigarFix() && !read->Alignments[scoreID].skip && read->Alignments[scoreID].cigarOpCount >= 0x10000);
	if (printLongCigar) { // write <readLength>S
		int clip_length = read->length;
		if (hardClip) clip_length = read->length - read->Alignments[scoreID].QStart - read->Alignments[scoreID].QEnd;
		Print("%dS\t", clip_length);
	} else {
		Print("%s\t", read->Alignments[scoreID].pBuffer1);
	}
	Print("%s\t", pRefName);//Ref. name of the mate/next fragment
	Print("%u\t", pLoc + report_offset);//Position of the mate/next fragment
	Print("%d\t", pDist);//observed Template LENgth

	if (hardClip) {
		Print("%.*s\t", read->length - read->Alignments[scoreID].QStart - read->Alignments[scoreID].QEnd,
				readseq + read->Alignments[scoreID].QStart);
	}
	else {
		Print("%.*s\t", read->length, readseq);
	}

	if (qltystr != 0) {
		if (hardClip)
		Print("%.*s\t", read->length - read->Alignments[scoreID].QStart - read->Alignments[scoreID].QEnd,
				qltystr + read->Alignments[scoreID].QStart);
		else
		Print("%.*s\t", read->length, qltystr);
	} else {
		Print("*\t");
	}

	//Optional fields
	if(rgId != 0) {
		Print("RG:Z:%s\t", rgId);
	}
	Print("AS:i:%d\t", (int) read->Scores[scoreID].Score.f);
	Print("NM:i:%d\t", read->Alignments[scoreID].NM);


	float identity = round(read->Alignments[scoreID].Identity * 10000.0f) / 10000.0f;
	Print("XI:f:%g\t", identity);
	Print("XS:i:%d\t", (int) 0.0f);
	Print("XE:i:%d\t", (int) read->Scores[scoreID].Score.f);
	Print("XR:i:%d\t", read->length - read->Alignments[scoreID].QStart - read->Alignments[scoreID].QEnd);
	Print("MD:Z:%s\t", read->Alignments[scoreID].pBuffer2);

	if(read->Alignments[scoreID].svType > -1) {
		Print("SV:i:%d\t", read->Alignments[scoreID].svType);
	}

	if(read->Calculated > 1) {
		bool first = true;
		for(int i = 0; i < read->Calculated; ++i) {
			if(i != scoreID && !read->Alignments[i].skip) {
				int refnamelen = 0;
				char const * refname = SequenceProvider.GetRefName(read->Scores[i].Location.getrefId(), refnamelen);

				char strand = '+';
				if (read->Scores[i].Location.isReverse()) {
					strand = '-';
				}
				if(first) {
					Print("SA:Z:");
					first = false;
				}
				Print("%.*s,%d,%c,%s,%d,%d;", refnamelen, refname, read->Scores[i].Location.m_Location + report_offset, strand, read->Alignments[i].pBuffer1, read->Alignments[i].MQ, read->Alignments[i].NM);
			}
		}
		if(!first) {
			Print("\t");
		}
	}

	Print("QS:i:%d\t", read->Alignments[scoreID].QStart);
	Print("QE:i:%d\t", read->length - read->Alignments[scoreID].QEnd);

	int clipped = read->Alignments[scoreID].QStart + read->Alignments[scoreID].QEnd;
	float covered = (read->length - clipped) * 100.0f / read->length;
	Print("CV:f:%f", covered);

	if (printLongCigar) { // write the real CIGAR at the CG:B,I tag
		Print("\tCG:B:I");
		char *p = read->Alignments[scoreID].pBuffer1;
		for (int i = 0; i < read->Alignments[scoreID].cigarOpCount; ++i) {
			long len = strtol(p, &p, 10);
			int op = 0;
			if (*p == 'M') op = 0;
			else if (*p == 'I') op = 1;
			else if (*p == 'D') op = 2;
			else if (*p == 'N') op = 3;
			else if (*p == 'S') op = 4;
			else if (*p == 'H') op = 5;
			else if (*p == '=') op = 7;
			else if (*p == 'X') op = 8;
			++p;
			unsigned int cigar1 = (unsigned int)len<<4 | op; // this is the binary representation of a CIGAR operation
			Print(",%d", cigar1);
		}
	}

	Print("\n");

}

void SAMWriter::DoWritePair(MappedRead const * const read1, int const scoreId1,
		MappedRead const * const read2, int const scoreId2) {
	//Proper pair
	int flags1 = 0x1;
	int flags2 = 0x1;
	if (read1->ReadId & 0x1) {
		//if read1 is the second pair
		flags1 |= 0x80;
		flags2 |= 0x40;
	} else {
		//if read1 is the first pair
		flags1 |= 0x40;
		flags2 |= 0x80;
	}
	if (!read1->hasCandidates() && !read2->hasCandidates()) {
		//Both mates unmapped
		DoWriteUnmappedRead(read2, flags2 | 0x8);
		DoWriteUnmappedRead(read1, flags1 | 0x8);
	} else if (!read1->hasCandidates()) {
		//First mate unmapped
		DoWriteReadGeneric(read2, scoreId2, "=",
				read2->Scores[scoreId2].Location.m_Location, 0,
				read2->mappingQlty, flags2 | 0x8);
		DoWriteUnmappedReadGeneric(read1,
				read2->Scores[scoreId2].Location.getrefId(), '=',
				read2->Scores[scoreId2].Location.m_Location,
				read2->Scores[scoreId2].Location.m_Location, 0, 0, flags1);
	} else if (!read2->hasCandidates()) {
		//Second mate unmapped
		DoWriteUnmappedReadGeneric(read2,
				read1->Scores[scoreId1].Location.getrefId(), '=',
				read1->Scores[scoreId1].Location.m_Location,
				read1->Scores[scoreId1].Location.m_Location, 0, 0, flags2);
		DoWriteReadGeneric(read1, scoreId1, "=",
				read1->Scores[scoreId1].Location.m_Location, 0,
				read1->mappingQlty, flags1 | 0x8);
	} else {
		if (!read1->HasFlag(NGMNames::PairedFail)) {
			//TODO: Check if correct!
			int distance = 0;

			flags1 |= 0x2;
			flags2 |= 0x2;
			if (!read1->Scores[scoreId1].Location.isReverse()) {
				distance = (read2->Scores[scoreId2].Location.m_Location
						+ read2->length - read2->Alignments[scoreId2].QStart
						- read2->Alignments[scoreId2].QEnd)
						- read1->Scores[scoreId1].Location.m_Location;
				DoWriteReadGeneric(read2, scoreId2, "=",
						read1->Scores[scoreId1].Location.m_Location,
						distance * -1, read2->mappingQlty, flags2);
				DoWriteReadGeneric(read1, scoreId1, "=",
						read2->Scores[scoreId2].Location.m_Location, distance,
						read1->mappingQlty, flags1 | 0x20);
			} else if (!read2->Scores[scoreId2].Location.isReverse()) {
				distance = (read1->Scores[scoreId1].Location.m_Location
						+ read1->length - read1->Alignments[scoreId1].QStart
						- read1->Alignments[scoreId1].QEnd)
						- read2->Scores[scoreId2].Location.m_Location;

				DoWriteReadGeneric(read2, scoreId2, "=",
						read1->Scores[scoreId1].Location.m_Location, distance,
						read2->mappingQlty, flags2 | 0x20);
				DoWriteReadGeneric(read1, scoreId1, "=",
						read2->Scores[scoreId2].Location.m_Location,
						distance * -1, read1->mappingQlty, flags1);
			}
		} else {
			int distance = 0;
			if (read1->Scores[scoreId1].Location.isReverse()) {
				flags2 |= 0x20;
			}
			if (read2->Scores[scoreId2].Location.isReverse()) {
				flags1 |= 0x20;
			}
			DoWriteReadGeneric(read2, scoreId2, SequenceProvider.GetRefName(read1->Scores[scoreId1].Location.getrefId(), distance), read1->Scores[scoreId1].Location.m_Location, 0, read2->mappingQlty, flags2);
			DoWriteReadGeneric(read1, scoreId1, SequenceProvider.GetRefName(read2->Scores[scoreId2].Location.getrefId(), distance), read2->Scores[scoreId2].Location.m_Location, 0, read1->mappingQlty, flags1);
		}
	}
	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer);
}

void SAMWriter::DoWriteUnmappedReadGeneric(MappedRead const * const read,
		int const refId, char const pRefName, int const loc, int const pLoc,
		int const pDist, int const mappingQlty, int flags = 0) {
	//SRR002320.10000027.1    4       *       0       0       *       *       0       0       TTTATGTTGTTAATGTGTTGGGTGAGTGCGCCCCAT    IIIIIIIIIIIIIIIIIIIIIIIIII

	NGM.AddUnmappedRead(read, 0);

	if(writeUnmapped) {
		NGM.AddWrittenRead(read->ReadId);
		flags |= 0x4;

		char const * readseq = read->Seq;
		int readlen = read->length;

		char * readname = read->name;

		char * qltystr = read->qlty;
		int qltylen = readlen;

		//mandatory fields
		Print("%s\t", readname);
		Print("%d\t", flags);
		if (refId > -1) {
			int refnamelen = 0;
			char const * refname = SequenceProvider.GetRefName(refId, refnamelen);
			Print("%.*s\t", refnamelen, refname);
		} else {
			Print("*\t");
		}
		Print("%d\t", loc + report_offset);

		Print("0\t");
		Print("*\t");
		Print("%c\t", pRefName); //Ref. name of the mate/next fragment
		Print("%d\t", pLoc + report_offset);//Position of the mate/next fragment
		Print("%d\t", pDist);//observed Template LENgth
		Print("%.*s\t", read->length, readseq);
		if (qltystr != 0) {
			Print("%.*s", qltylen, qltystr);
		} else {
			Print("*");
		}

		if(rgId != 0) {
			Print("\tRG:Z:%s", rgId);
		}

		Print("\n");
	}
}

void SAMWriter::DoWriteUnmappedRead(MappedRead const * const read, int flags) {
	DoWriteUnmappedReadGeneric(read, -1, '*', -1, -1, 0, 0, flags | 0x04);
	m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer);
}

void SAMWriter::DoWriteEpilog() {

}
