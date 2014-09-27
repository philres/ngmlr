#include "BAMWriter.h"

#ifdef _BAM

#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "Config.h"
#include "SequenceProvider.h"
#include "NGM.h"

#include "Version.h"

using namespace BamTools;

void BAMWriter::DoWriteProlog() {
	//TODO: check correct format	;

	//NGMLock(&m_OutputMutex);
	SamHeader header;
	RefVector refs;

	SamProgram program;

	header.Version = "1.0";
	header.SortOrder = "unsorted";

	program.ID = "ngm";
	program.Name = "ngm";
	std::stringstream version;
	version << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_BUILD;
	program.Version = version.str();
	std::stringstream cmdline;
	cmdline << "\"" << std::string(Config.GetString("cmdline")) << "\"";
	program.CommandLine = cmdline.str();

	header.Programs.Add(program);

	if (RG.size() > 0) {
		SamReadGroup rg(RG);

		if (Config.Exists(RG_CN))
			rg.SequencingCenter = std::string(Config.GetString(RG_CN));

		if (Config.Exists(RG_DS))
			rg.Description = std::string(Config.GetString(RG_DS));

		if (Config.Exists(RG_DT))
			rg.ProductionDate = std::string(Config.GetString(RG_DT));

		if (Config.Exists(RG_FO))
			rg.FlowOrder = std::string(Config.GetString(RG_FO));

		if (Config.Exists(RG_KS))
			rg.KeySequence = std::string(Config.GetString(RG_KS));

		if (Config.Exists(RG_LB))
			rg.Library = std::string(Config.GetString(RG_LB));

		if (Config.Exists(RG_PG))
			rg.Program = std::string(Config.GetString(RG_PG));

		if (Config.Exists(RG_PI))
			rg.PredictedInsertSize = std::string(Config.GetString(RG_PI));

		if (Config.Exists(RG_PL))
			rg.SequencingTechnology = std::string(Config.GetString(RG_PL));

		if (Config.Exists(RG_PU))
			rg.PlatformUnit = std::string(Config.GetString(RG_PU));

		if (Config.Exists(RG_SM))
			rg.Sample = std::string(Config.GetString(RG_SM));

		header.ReadGroups.Add(rg);
	}

	char const * refName = 0;
	int refNameLength = 0;
	int ref = Config.GetInt("ref_mode");

	if (ref == -1)
		for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
			refName = SequenceProvider.GetRefName(i, refNameLength);
			//Print("@SQ\tSN:%.*s\tLN:%d\n", refNameLength, refName, SequenceProvider.GetRefLen(i));
			RefData bamRef(std::string(refName, refNameLength),SequenceProvider.GetRefLen(i));
			refs.push_back(bamRef);
			if (NGM.DualStrand())
			++i;
		}
		else {
			refName = SequenceProvider.GetRefName(ref * ((NGM.DualStrand()) ? 2 : 1), refNameLength);
			RefData bamRef(std::string(refName, refNameLength), SequenceProvider.GetRefLen(ref * ((NGM.DualStrand()) ? 2 : 1)));
			refs.push_back(bamRef);
		}

	if (!writer->Open(file, header, refs)) {
		Log.Error("Could not open output BAM file");
		return;
	}
	//NGMUnlock(&m_OutputMutex);
}

void BAMWriter::translate_flag(BamAlignment * al, int flags) {
	al->SetIsPaired(flags & 0x1);
	al->SetIsProperPair(flags & 0x2);
	al->SetIsMapped(!(flags & 0x4));
	al->SetIsMateMapped(!(flags & 0x8));
	al->SetIsMateReverseStrand(flags & 0x20);
	al->SetIsFirstMate(flags & 0x40);
	al->SetIsSecondMate(flags & 0x80);
}

void BAMWriter::addAdditionalInfo(const MappedRead * const read, BamAlignment * al) {

	std::stringstream ss(read->AdditionalInfo);
	std::string token;
	char delim = '\t';
	while (std::getline(ss, token, delim)) {
		if (token.length() > 5) {
			int first = token.find(':');
			int last = token.rfind(':');
			if (first > 0 && last > 0 && first < last) {
				if (token.substr(first + 1, last - first - 1).size() == 1) {
					char type = token.substr(first + 1, last - first - 1)[0];
					if (type == Constants::BAM_TAG_TYPE_INT8 || type == Constants::BAM_TAG_TYPE_INT16
							|| type == Constants::BAM_TAG_TYPE_INT32) {
						al->AddTag(token.substr(0, first), token.substr(first + 1, last - first - 1),
								atoi(token.substr(last + 1, token.length() - last).c_str()));
					} else if (type == Constants::BAM_TAG_TYPE_UINT8 || type == Constants::BAM_TAG_TYPE_UINT16
							|| type == Constants::BAM_TAG_TYPE_UINT32) {
						token.substr(first + 1, last - first - 1), atoi(token.substr(last + 1, token.length() - last).c_str());
					} else if (type == Constants::BAM_TAG_TYPE_STRING || type == Constants::BAM_TAG_TYPE_ASCII) {
						al->AddTag(token.substr(0, first), token.substr(first + 1, last - first - 1),
								token.substr(last + 1, token.length() - last));
					} else if (type == Constants::BAM_TAG_TYPE_FLOAT) {
						token.substr(first + 1, last - first - 1), atof(token.substr(last + 1, token.length() - last).c_str());
					}
				}
			}
		}
	}
}

void BAMWriter::DoWriteReadGeneric(MappedRead const * const read, int const scoreId, int const pRef, int const pLoc, int const pDist, int const mappingQlty, int flags) {
	NGM.AddWrittenRead(read->ReadId);
	static bool const hardClip = Config.GetInt("hard_clip", 0, 1) == 1 || Config.GetInt("silent_clip", 0, 1) == 1;

	BamAlignment * al = new BamAlignment();

	char const * readseq = read->Seq;
	int readlen = read->length;
	char * readname = read->name;

	char * qltystr = read->qlty;
	int qltylen = readlen;

	al->AlignmentFlag = 0;
	translate_flag(al, flags);
	al->SetIsMapped(true);

	if (read->Scores[scoreId].Location.isReverse()) {
		readseq = read->RevSeq;
		if (qltystr != 0 && qltystr[0] != '*') {
			std::reverse(qltystr, &qltystr[read->length]);
		}
		al->SetIsReverseStrand(true);
	}

//al->Name = std::string(readname, readnamelen);
	al->Name = std::string(readname);
	al->Length = readlen;
	al->MapQuality = read->mappingQlty;
	al->Position = read->Scores[scoreId].Location.m_Location;

	if (hardClip) {
		al->QueryBases = std::string(readseq + read->Alignments[scoreId].QStart,
				read->length - read->Alignments[scoreId].QStart - read->Alignments[scoreId].QEnd);
	}
	else {
		al->QueryBases = std::string(readseq, readlen);
	}

	al->RefID = read->Scores[scoreId].Location.getrefId() / 2;

//Paired
//if (pRefName == '=')
//al->MateRefID = read->TLS()->Location.m_RefId / 2;
	if (pRef > -1)
	al->MateRefID = pRef / 2;
	al->MatePosition = pLoc;
	al->InsertSize = pDist;

	int length = atoi(read->Alignments[scoreId].pBuffer1);
	int i = 0;
//Log.Message("%s: %s", read->name, read->Buffer1);
	while (read->Alignments[scoreId].pBuffer1[i] != '\0') {
		if (read->Alignments[scoreId].pBuffer1[i] < '0' || read->Alignments[scoreId].pBuffer1[i] > '9') {
			al->CigarData.push_back(CigarOp(read->Alignments[scoreId].pBuffer1[i], length));
			length = atoi(read->Alignments[scoreId].pBuffer1 + i + 1);
		}
		i += 1;
	}

	if (qltystr != 0 && qltylen > 0 && qltylen == readlen) {
		if (hardClip)
		al->Qualities = std::string(qltystr + read->Alignments[scoreId].QStart,
				read->length - read->Alignments[scoreId].QStart - read->Alignments[scoreId].QEnd);
		else
		al->Qualities = std::string(qltystr, qltylen);
	} else {
		if (hardClip)
		al->Qualities = std::string(read->length - read->Alignments[scoreId].QStart - read->Alignments[scoreId].QEnd, ':');
		else
		al->Qualities = std::string(readlen, ':');
	}

	if (read->AdditionalInfo != 0) {
		addAdditionalInfo(read, al);
	}

//	//Optional fields
	al->AddTag("AS", "i", (int) read->Scores[scoreId].Score.f);
	al->AddTag("NM", "i", read->Alignments[scoreId].NM);

	if (Config.GetInt("bs_mapping") == 1) {
		if (!(read->ReadId & 1) || read->Paired == 0) {
			if (read->Scores[scoreId].Location.isReverse()) {
				al->AddTag("ZS", "Z", std::string("-+"));
			} else {
				al->AddTag("ZS", "Z", std::string("++"));
			}
		} else {
			if (read->Scores[scoreId].Location.isReverse()) {
				al->AddTag("ZS", "Z", std::string("+-"));
			} else {
				al->AddTag("ZS", "Z", std::string("--"));
			}
		}
	}

	float identity = round(read->Alignments[scoreId].Identity * 10000.0f) / 10000.0f;
	al->AddTag("XI", "f", identity);
	al->AddTag("X0", "i", (int) read->numTopScores);
//TODO: fix. Calculated used to be the number of score computed. Now it is the number of computed alignments.
//Thus it can't be used for X1 anymore.
//al->AddTag("X1", "i", (int) (read->Calculated - read->EqualScoringCount));
	al->AddTag("XE", "i", (int) read->s);
	al->AddTag("XR", "i", read->length - read->Alignments[scoreId].QStart - read->Alignments[scoreId].QEnd);
	al->AddTag("MD", "Z", std::string(read->Alignments[scoreId].pBuffer2));

	if(RG.size() > 0) {
		al->AddTag("RG", "Z", RG);
	}

	buffer[bufferIndex++] = al;
	if (bufferIndex == (10000 - 1)) {
		writer->SaveAlignment(buffer, bufferIndex);
		bufferIndex = 0;
	}

//NGMUnlock(&m_OutputMutex);
}

void BAMWriter::DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, int const pRef, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags) {

	NGM.AddUnmappedRead(read, MFAIL_NOCAND);

	if(writeUnmapped) {
		NGM.AddWrittenRead(read->ReadId);
		BamAlignment * al = new BamAlignment();

		char const * readseq = read->Seq;
		int readlen = read->length;

		char * readname = read->name;
		//int readnamelen = read->nameLength;

		char * qltystr = read->qlty;
		int qltylen = readlen;

		flags |= 0x4;
		al->AlignmentFlag = 0;
		translate_flag(al, flags);

		//al->Name = std::string(readname, readnamelen);
		al->Name = std::string(readname);
		al->Length = readlen;
		al->MapQuality = 0;

		if(loc > -1) {
			al->Position = loc;
		}

		if(refId > -1) {
			al->RefID = refId / 2;
		}

		if (pRef > -1) {
			al->MateRefID = pRef / 2;
		}
		if(pLoc > -1) {
			al->MatePosition = pLoc;
		}

		std::string seq(readseq, readlen);
		al->QueryBases = seq;

		if (qltystr != 0 && qltylen > 0 && qltylen == readlen) {
			al->Qualities = std::string(qltystr, qltylen);
		} else {
			al->Qualities = std::string(readlen, ':');
		}

		if(RG.size() > 0) {
			al->AddTag("RG", "Z", RG);
		}

		//	//Optional fields
//	al->AddTag("AS", "i", (int) read->TLS()->Score.f);
//	al->AddTag("MD", "Z", std::string(read->Buffer2));
//	al->AddTag("X0", "i", (int) read->EqualScoringCount);
//	al->AddTag("XI", "f", read->Identity);

//writer->SaveAlignment(al);
		buffer[bufferIndex++] = al;
		if (bufferIndex == (10000 - 1)) {
			writer->SaveAlignment(buffer, bufferIndex);
			bufferIndex = 0;
		}

		if (read->AdditionalInfo != 0) {
			addAdditionalInfo(read, al);
		}
	}
}

void BAMWriter::DoWriteRead(MappedRead const * const read, int const scoreId) {
	DoWriteReadGeneric(read, scoreId, -1, -1, 0, read->mappingQlty, 0);
}

void BAMWriter::DoWriteUnmappedRead(MappedRead const * const read, int flags) {
	DoWriteUnmappedReadGeneric(read, -1, -1, -1, -1, 0, 0, flags);
}

void BAMWriter::DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2) {
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
		DoWriteReadGeneric(read2, scoreId2, read2->Scores[scoreId2].Location.getrefId(), read2->Scores[scoreId2].Location.m_Location, 0,
				read2->mappingQlty, flags2 | 0x8);
		DoWriteUnmappedReadGeneric(read1, read2->Scores[scoreId2].Location.getrefId(), read2->Scores[scoreId2].Location.getrefId(),
				read2->Scores[scoreId2].Location.m_Location, read2->Scores[scoreId2].Location.m_Location, 0, 0, flags1);
	} else if (!read2->hasCandidates()) {
		//Second mate unmapped
		DoWriteUnmappedReadGeneric(read2, read1->Scores[scoreId1].Location.getrefId(), read1->Scores[scoreId1].Location.getrefId(),
				read1->Scores[scoreId1].Location.m_Location, read1->Scores[scoreId1].Location.m_Location, 0, 0, flags2);
		DoWriteReadGeneric(read1, scoreId1, read1->Scores[scoreId1].Location.getrefId(), read1->Scores[scoreId1].Location.m_Location, 0,
				read1->mappingQlty, flags1 | 0x8);
	} else {
		if (!read1->HasFlag(NGMNames::PairedFail)) {
			//TODO: Check if correct!
			int distance = 0;
			flags1 |= 0x2;
			flags2 |= 0x2;
			if (!read1->Scores[scoreId1].Location.isReverse()) {
				distance = read2->Scores[scoreId2].Location.m_Location + read2->length - read1->Scores[scoreId1].Location.m_Location;
				DoWriteReadGeneric(read2, scoreId2, read2->Scores[scoreId2].Location.getrefId(),
						read1->Scores[scoreId1].Location.m_Location, distance * -1, read2->mappingQlty, flags2);
				DoWriteReadGeneric(read1, scoreId1, read2->Scores[scoreId2].Location.getrefId(),
						read2->Scores[scoreId2].Location.m_Location, distance, read1->mappingQlty, flags1 | 0x20);
			} else if (!read2->Scores[scoreId2].Location.isReverse()) {
				distance = read1->Scores[scoreId1].Location.m_Location + read1->length - read2->Scores[scoreId2].Location.m_Location;

				DoWriteReadGeneric(read2, scoreId2, read2->Scores[scoreId2].Location.getrefId(),
						read1->Scores[scoreId1].Location.m_Location, distance, read2->mappingQlty, flags2 | 0x20);
				DoWriteReadGeneric(read1, scoreId1, read2->Scores[scoreId2].Location.getrefId(),
						read2->Scores[scoreId2].Location.m_Location, distance * -1, read1->mappingQlty, flags1);
			}
		} else {
			if (read1->Scores[scoreId1].Location.isReverse()) {
				flags2 |= 0x20;
			}
			if (read2->Scores[scoreId2].Location.isReverse()) {
				flags1 |= 0x20;
			}
			DoWriteReadGeneric(read2, scoreId2, read1->Scores[scoreId1].Location.getrefId(), read1->Scores[scoreId1].Location.m_Location, 0,
					read2->mappingQlty, flags2);
			DoWriteReadGeneric(read1, scoreId1, read2->Scores[scoreId2].Location.getrefId(), read2->Scores[scoreId2].Location.m_Location, 0,
					read1->mappingQlty, flags1);
		}
	}
}
void BAMWriter::DoWriteEpilog() {
//	writer->Close();
}

#endif
