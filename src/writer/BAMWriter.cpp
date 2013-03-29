#include "BAMWriter.h"

#ifdef _BAM

#include <string.h>

#include "api/BamWriter.h"
#include "api/SamHeader.h"

#include "Config.h"
#include "SequenceProvider.h"
#include "NGM.h"

using namespace BamTools;

//using BamTools::BamWriter;
//using BamTools::SamHeader;
//Format: http://samtools.sourceforge.net/SAM1.pdf

//static inline void rev(char * s);

BamWriter writer;

void BAMWriter::DoWriteProlog() {
	//TODO: check correct format	;

	NGMLock(&m_OutputMutex);
	SamHeader header;
	RefVector refs;

	SamProgram program;

	header.Version = "1.0";
	header.SortOrder = "unsorted";

	program.Name = "ngm";
	program.Version = "0.0.1";
	program.CommandLine = std::string(Config.GetString("cmdline"));

	header.Programs.Add(program);

	char const * refName = 0;
	int refNameLength = 0;
	int ref = Config.GetInt("ref_mode");

	if (ref == -1)
	for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
		refName = SequenceProvider.GetRefName(i, refNameLength);
		//Print("@SQ\tSN:%.*s\tLN:%d\n", refNameLength, refName, SequenceProvider.GetRefLen(i));
		RefData bamRef(std::string(refName, refNameLength), SequenceProvider.GetRefLen(i));
		refs.push_back(bamRef);
		if (NGM.DualStrand())
		++i;
	}
	else {
		refName = SequenceProvider.GetRefName(ref * ((NGM.DualStrand()) ? 2 : 1), refNameLength);
		RefData bamRef(std::string(refName, refNameLength), SequenceProvider.GetRefLen(ref * ((NGM.DualStrand()) ? 2 : 1)));
		refs.push_back(bamRef);
	}

	if (!writer.Open(file, header, refs)) {
		Log.Error("Could not open output BAM file");
		return;
	}
	NGMUnlock(&m_OutputMutex);
}

void BAMWriter::translate_flag(BamAlignment &al, int flags) {
	al.SetIsPaired(flags & 0x1);
	al.SetIsProperPair(flags & 0x2);
	al.SetIsMateMapped(!(flags & 0x8));
	al.SetIsMateReverseStrand(flags & 0x20);
	al.SetIsFirstMate(flags & 0x40);
	al.SetIsSecondMate(flags & 0x80);
}

void BAMWriter::DoWriteReadGeneric(MappedRead const * const read, int const pRef, int const pLoc, int const pDist,
		int const mappingQlty, int flags) {

	NGMLock(&m_OutputMutex);
	static bool const hardClip = Config.GetInt("hard_clip", 0, 1) == 1 || Config.GetInt("silent_clip", 0, 1) == 1;

	BamAlignment al;

	char const * readseq = read->Seq;
	int readlen = read->length;

	char * readname = read->name;
//	int readnamelen = read->nameLength;

	char * qltystr = read->qlty;
	int qltylen = readlen;

	al.AlignmentFlag = 0;
	translate_flag(al, flags);
	al.SetIsMapped(true);

	if ((read->Strand == '-')) {
		readseq = read->RevSeq;
		al.SetIsReverseStrand(true);
	}

	//al.Name = std::string(readname, readnamelen);
	al.Name = std::string(readname);
	al.Length = readlen;
	al.MapQuality = read->mappingQlty;
	al.Position = read->TLS()->Location.m_Location;

	if (hardClip)
	al.QueryBases = std::string(readseq + read->QStart, read->length - read->QStart - read->QEnd);
	else
	al.QueryBases = std::string(readseq, readlen);

	al.RefID = read->TLS()->Location.m_RefId / 2;

	//Paired
	//if (pRefName == '=')
	//al.MateRefID = read->TLS()->Location.m_RefId / 2;
	if(pRef > -1)
		al.MateRefID = pRef / 2;
	al.MatePosition = pLoc;
	al.InsertSize = pDist;

	int length = atoi(read->Buffer1);
	int i = 0;
	//Log.Message("%s: %s", read->name, read->Buffer1);
	while (read->Buffer1[i] != '\0') {
		if (read->Buffer1[i] < '0' || read->Buffer1[i] > '9') {
			al.CigarData.push_back(CigarOp(read->Buffer1[i], length));
			length = atoi(read->Buffer1 + i + 1);
		}
		i += 1;
	}

	if (qltystr != 0 && qltylen > 0 && qltylen == readlen) {
		if (hardClip)
		al.Qualities = std::string(qltystr + read->QStart, read->length - read->QStart - read->QEnd);
		else
		al.Qualities = std::string(qltystr, qltylen);
	} else {
		if (hardClip)
		al.Qualities = std::string(read->length - read->QStart - read->QEnd, ':');
		else
		al.Qualities = std::string(readlen, ':');
	}

	//	//Optional fields
	al.AddTag("AS", "i", (int) read->TLS()->Score.f);
	al.AddTag("NM", "i", read->NM);

	if (Config.GetInt("bs_mapping") == 1) {
		if (!(read->ReadId & 1)) {
			if (read->Strand == '-') {
				al.AddTag("ZS", "Z", std::string("-+"));
			} else {
				al.AddTag("ZS", "Z", std::string("++"));
			}
		} else {
			if (read->Strand == '-') {
				al.AddTag("ZS", "Z", std::string("+-"));
			} else {
				al.AddTag("ZS", "Z", std::string("--"));
			}
		}
	}

	al.AddTag("XI", "f", read->Identity);
	al.AddTag("X0", "i", (int) read->EqualScoringCount);
	al.AddTag("X1", "i", (int) (read->Calculated - read->EqualScoringCount));
	al.AddTag("XE", "i", (int) read->s);
	al.AddTag("XR", "i", read->length - read->QStart - read->QEnd);
	al.AddTag("MD", "Z", std::string(read->Buffer2));

	writer.SaveAlignment(al);
	NGMUnlock(&m_OutputMutex);
}

void BAMWriter::DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, char const pRefName, int const loc,
		int const pLoc, int const pDist, int const mappingQlty, int flags) {

	NGMLock(&m_OutputMutex);
	NGM.AddUnmappedRead(read, MFAIL_NOCAND);

	BamAlignment al;

	char const * readseq = read->Seq;
	int readlen = read->length;

	char * readname = read->name;
	//int readnamelen = read->nameLength;

	char * qltystr = read->qlty;
	int qltylen = readlen;

	al.AlignmentFlag = 0;
	translate_flag(al, flags);
	al.SetIsMapped(false);

	//al.Name = std::string(readname, readnamelen);
	al.Name = std::string(readname);
	al.Length = readlen;
	al.MapQuality = 0;

	std::string seq(readseq, readlen);
	al.QueryBases = seq;

	//al.Position = 0;
	//al.RefID = 0;

	//Log.Message("%s: %s", read->name, read->Buffer1);

	if (qltystr != 0 && qltylen > 0 && qltylen == readlen) {
		al.Qualities = std::string(qltystr, qltylen);
	} else {
		al.Qualities = std::string(readlen, ':');
	}

	//	//Optional fields
//	al.AddTag("AS", "i", (int) read->TLS()->Score.f);
//	al.AddTag("MD", "Z", std::string(read->Buffer2));
//	al.AddTag("X0", "i", (int) read->EqualScoringCount);
//	al.AddTag("XI", "f", read->Identity);

	writer.SaveAlignment(al);
	NGMUnlock(&m_OutputMutex);
}

void BAMWriter::DoWriteRead(MappedRead const * const read) {
	DoWriteReadGeneric(read, -1, -1, 0, read->mappingQlty);
}

void BAMWriter::DoWriteUnmappedRead(MappedRead const * const read, int flags) {
	DoWriteUnmappedReadGeneric(read, -1, '*', -1, -1, 0, 0, flags);
}

void BAMWriter::DoWritePair(MappedRead const * const read1, MappedRead const * const read2) {
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
		DoWriteUnmappedRead(read2, flags2);
		DoWriteUnmappedRead(read1, flags1);
	} else if (!read1->hasCandidates()) {
		//First mate unmapped
		flags2 |= 0x8;

		DoWriteReadGeneric(read2, -1, read2->TLS()->Location.m_Location, 0, read2->mappingQlty, flags2);
		DoWriteUnmappedReadGeneric(read1, read2->TLS()->Location.m_RefId, read2->TLS()->Location.m_RefId, read2->TLS()->Location.m_Location,
				read2->TLS()->Location.m_Location, 0, 0, flags1);
	} else if (!read2->hasCandidates()) {
		flags1 |= 0x8;
		//Second mate unmapped
		DoWriteUnmappedReadGeneric(read2, read1->TLS()->Location.m_RefId, read1->TLS()->Location.m_RefId, read1->TLS()->Location.m_Location,
				read1->TLS()->Location.m_Location, 0, 0, flags2);
		DoWriteReadGeneric(read1, -1, read1->TLS()->Location.m_Location, 0, read1->mappingQlty, flags1);
	} else {
		if (!read1->HasFlag(NGMNames::PairedFail)) {
			//TODO: Check if correct!
			int distance = 0;
			flags1 |= 0x2;
			flags2 |= 0x2;
			if (read1->Strand == '+') {
				distance = read2->TLS()->Location.m_Location + read2->length - read1->TLS()->Location.m_Location;
				DoWriteReadGeneric(read2, read2->TLS()->Location.m_RefId, read1->TLS()->Location.m_Location, distance * -1, read2->mappingQlty, flags2);
				DoWriteReadGeneric(read1, read2->TLS()->Location.m_RefId, read2->TLS()->Location.m_Location, distance, read1->mappingQlty, flags1 | 0x20);
			} else if (read2->Strand == '+') {
				distance = read1->TLS()->Location.m_Location + read1->length - read2->TLS()->Location.m_Location;

				DoWriteReadGeneric(read2, read2->TLS()->Location.m_RefId, read1->TLS()->Location.m_Location, distance, read2->mappingQlty, flags2 | 0x20);
				DoWriteReadGeneric(read1, read2->TLS()->Location.m_RefId, read2->TLS()->Location.m_Location, distance * -1, read1->mappingQlty, flags1);
			}
		} else {
			if (read1->Strand == '-') {
				flags2 |= 0x20;
			}
			if (read2->Strand == '-') {
				flags1 |= 0x20;
			}
			DoWriteReadGeneric(read2, read1->TLS()->Location.m_RefId,
					read1->TLS()->Location.m_Location, 0, read2->mappingQlty, flags2);
			DoWriteReadGeneric(read1, read2->TLS()->Location.m_RefId,
					read2->TLS()->Location.m_Location, 0, read1->mappingQlty, flags1);

			//DoWriteRead(read2);
			//DoWriteRead(read1);
			//DoWriteUnmappedRead(read2);
			//DoWriteUnmappedRead(read1);
		}
	}
}
void BAMWriter::DoWriteEpilog() {
	writer.Close();
}

#endif
