#include "FFormatWriter.h"

#include "Config.h"
#include "SequenceProvider.h"
#include "NGM.h"

#include <string>
std::string add_timestamp(std::string);

void FFormatWriter::DoWriteProlog() {
	char const * refName = 0;
	int len;
	Print("%s", add_timestamp("NGM Output (generated %s)\n\n").c_str());

	int ref = Config.GetInt("ref_mode");
	int cnt = SequenceProvider.GetRefCount() >> (NGM.DualStrand() ? 1 : 0);
	Print("Reference sequence/chromosome count (used for search/total in file): %i / %i\n", (ref == -1) ? cnt : 1, cnt);

	if (ref == -1)
		for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
			refName = SequenceProvider.GetRefName(i, len);
			Print("> %.*s\n", len, refName);
			Print("  Length: %i\n", SequenceProvider.GetRefLen(i));

			if (NGM.DualStrand())
				++i;
		}
	else {
		refName = SequenceProvider.GetRefName(ref * ((NGM.DualStrand()) ? 2 : 1), len);
		Print("> %.*s\n", len, refName);
	}
	Print("\n");
}

void FFormatWriter::DoWriteRead(MappedRead const * const read) {
	static int const report_offset = 0;

	char * readname = read->name;
	uint readnamelen = read->nameLength;

	int refnamelen = 0;
	char const * refname = SequenceProvider.GetRefName(read->TLS()->Location.m_RefId, refnamelen);

	int refnamelenmod = refnamelen;
//	if (refnamelenmod >= 4) // Wenn refname mehr als 3 Zeichen lang ist, wird beim ersten Leerzeichen abgebrochen
//	{
//		refnamelenmod = 4;
//		while (refnamelenmod < refnamelen && *(refname+refnamelenmod) != ' ')
//			++refnamelenmod;
//	}

//char * qltystr = read->qlty;
//int qltylen = read->length;
//SequenceProvider.GetQualityString(qltystr, read->ReadId, qltylen);

	Print("ID: %i.%i\tRead: %.*s\tDir: %c%c\tPos: %i\tChr: %.*s\tScore: %.1f\n", read->ReadId, read->EqualScoringID, readnamelen, readname,
			read->Strand, read->Strand, read->TLS()->Location.m_Location + report_offset, refnamelenmod, refname, read->TLS()->Score.f);
	Print("times: %i\tidentity: %f\tQ_start: %i\tQ_end: %i\n", read->EqualScoringCount, read->Identity, read->QStart, read->QEnd);

	// �berspringe Leerzeichen am Anfang
	int space_skip = 0;
	while (read->Buffer1[space_skip] == ' ')
		++space_skip;

	Print("%s\n%s\n", read->Buffer1 + space_skip, read->Buffer2 + space_skip);

//	// Output Quality-String (+gaps)
//	if (qltylen > 0)
//	{
//		// TODO: QStart, QEnd, Strand reincodieren
//		int pos = 0;
//		int qpos = 0;
//		while (read->Buffer2[space_skip+pos] != 0)
//		{
//			if (read->Buffer2[space_skip+pos] == '-')
//				Print("-");
//			else
//				Print("%c", *(qltystr+qpos++));
//			++pos;
//		}
//	}

	Print("\n\n");
}

void FFormatWriter::DoWritePair(MappedRead const * const read1, MappedRead const * const read2) {
	if (read1->hasCandidates()) {
		DoWriteRead(read1);
	} else {
		DoWriteUnmappedRead(read1);
	}
	if (read2->hasCandidates()) {
		DoWriteRead(read2);
	} else {
		DoWriteUnmappedRead(read2);
	}
}

void FFormatWriter::DoWriteUnmappedRead(MappedRead const * const read) {

}

void FFormatWriter::DoWriteEpilog() {
}

