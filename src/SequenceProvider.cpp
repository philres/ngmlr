#include "SequenceProvider.h"
#include "NGM.h"

#include <vector>
#include <map>
#include <string>
#include <cmath>

#include <string.h>
#include <limits.h>

#include "Debug.h"

#include <iostream>

#undef module_name
#define module_name "SEQPROV"

_SequenceProvider * _SequenceProvider::pInstance = 0;

_SequenceProvider & _SequenceProvider::Instance() {
	if (pInstance == 0)
		pInstance = new _SequenceProvider();

	return *pInstance;
}

//static int const maxRefCount = 32768;
static int const maxRefCount = 32768000;
static uint const refEncCookie = 0x74656;

//const size_t _SequenceProvider::MAX_REF_NAME_LENGTH = 150;

//bool checkChar(char c) {
//	return (c == 'A') | (c == 'C') | (c == 'T') | (c == 'G');
//}
//
//// checked die sequence seq auf ACGT, liefert das erste auftreten eines von diesen
//// buchstaben unterschiedlichen zeichens ...4-basen-worte (32bit) verwenden?
//int checkSequence(char * seq, uint len) {
//	for (uint n = 0; n < len; ++n) {
//		if (!checkChar(*(seq + n)))
//			return n;
//		++n;
//	}
//	return -1;
//}

std::string CheckFile(std::string filename, char const * const name) {
	if (!FileExists(filename.c_str())) {
		Log.Error("%s file not found (%s)", name, filename.c_str());
		Fatal();
	}
	return filename;
}

#include <stdio.h>
#include <zlib.h>

#include "Timing.h"

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//static inline char enc4(char c) {
//	c = toupper(c);
//	switch (c) {
//	case 'A':
//		return 0x1;
//	case 'T':
//		return 0x4;
//	case 'G':
//		return 0x7;
//	case 'C':
//		return 0x3;
//	}
//	return 0xE;
//}
//
//static inline char dec4High(unsigned char c) {
//	c = (c >> 4) | 0x40;
//	if (c == 0x44)
//		c |= 0x10;
//	return c;
//}
//
//static inline char dec4Low(unsigned char c) {
//	c = (c & 0xF) | 0x40;
//	if (c == 0x44)
//		c |= 0x10;
//	return c;
//}

static inline char enc4(char c) {
	c = toupper(c);
	switch (c) {
		case 'A':
			return 0;
		case 'T':
			return 1;
		case 'G':
			return 2;
		case 'C':
			return 3;
	}
	return 4;
}

static inline char dec4(char c) {
	switch (c) {
		case 0:
			return 'A';
		case 1:
			return 'T';
		case 2:
			return 'G';
		case 3:
			return 'C';
		case 4:
			return 'N';
	}
	throw "Error in ref encoding!";
}

static inline char dec4High(unsigned char c) {
	return dec4(c >> 4);
}

static inline char dec4Low(unsigned char c) {
	return dec4(c & 0xF);
}

SequenceLocation _SequenceProvider::convert(MappedRead * read, uint m_Location) {
	SequenceLocation loc;
	int refCount = SequenceProvider.GetRefCount();

	int j = 0;
	while (j < refCount && m_Location >= SequenceProvider.GetRefStart(j)) {
		j += (NGM.DualStrand()) ? 2 : 1;
	}
	if (j == 0) {
		//Log.Message("%s (%d) - %s: %s, %s", read->name, read->ReadId, read->Seq, read->Buffer1, read->Buffer2);
		Log.Error("Couldn't resolve mapping position: %u!", m_Location);
		Fatal();
	}
	j -= (NGM.DualStrand()) ? 2 : 1;

	loc.m_Location = m_Location - SequenceProvider.GetRefStart(j);
	loc.setRefId(j);

	return loc;
}

int _SequenceProvider::readEncRefFromFile(char const * fileName) {
	Log.Message("Reading encoded reference from %s", fileName);
	Timer wtmr;
	wtmr.ST();

	uint encRefSize = 0;
	uint refCount = 0;
	uint cookie = 0;

	FILE *fp;
	fp = fopen(fileName, "rb");

	fread(&cookie, sizeof(uint), 1, fp);
	fread(&refCount, sizeof(uint), 1, fp);
	fread(&binRefIndex, sizeof(uint), 1, fp);
	fread(&encRefSize, sizeof(uint), 1, fp);
	if (cookie != refEncCookie) {
		fclose(fp);
		Log.Error("Invalid encoded reference file found: %s.", fileName);
		Log.Error("Please delete it and run NGM again.");
		Fatal();
	}
	if(refCount > maxRefCount) {
		Log.Error("Currently NextGenMap can't handle more than %d reference sequences.", maxRefCount);
		Fatal();
	}
	binRefIdx = new RefIdx[refCount];
	fread(binRefIdx, sizeof(RefIdx), refCount, fp);

	binRef = new char[encRefSize];
	fread(binRef, sizeof(char), encRefSize, fp);
	fclose(fp);
	Log.Message("Reading from disk took %.2fs", wtmr.ET());

	return refCount;
}

void _SequenceProvider::writeEncRefToFile(char const * fileName, uint const refCount, uint const encRefSize) {
	if (!Config.GetInt("skip_save")) {
		Timer wtmr;
		wtmr.ST();
		Log.Message("Writing encoded reference to %s", fileName);
		FILE *fp;
		fp = fopen(fileName, "wb");
		fwrite(&refEncCookie, sizeof(uint), 1, fp);
		fwrite(&refCount, sizeof(uint), 1, fp);
		fwrite(&binRefIndex, sizeof(uint), 1, fp);
		fwrite(&encRefSize, sizeof(uint), 1, fp);

		fwrite(binRefIdx, sizeof(RefIdx), refCount, fp);
		fwrite(binRef, sizeof(char), encRefSize, fp);

		fclose(fp);
		Log.Message("Writing to disk took %.2fs", wtmr.ET());
	}
}

long getSize(char const * const file) {
	gzFile gzfp;
	kseq_t *seq;
	gzfp = gzopen(Config.GetString("ref"), "r");
	seq = kseq_init(gzfp);
	int l = 0;
	//1000 -> padding at beginning
	long size = 1000;
	while ((l = kseq_read(seq)) >= 0) {
		//1000 -> 1000 x N padding
		int s = (seq->seq.l | 1) + 1;
		size += s + 1000;
	}
	kseq_destroy(seq);
	gzclose(gzfp);
	return size;
}

void _SequenceProvider::Init() {
	Log.Message("Init sequence provider.");
	if (!Config.Exists("ref")) {
		Log.Error("No reference file specified.");
		Fatal();
	}
//	if (NGM.Paired() && !Config.Exists("pqry")) {
//		Log.Error("No paired query file specified.");
//		Fatal();
//	}

	CheckFile(refBaseFileName = Config.GetString("ref"), "RefBase");

	refFileName = Config.GetString("ref") + std::string("-enc.ngm");

	if (FileExists(refFileName.c_str())) {
		//Read
		refCount = readEncRefFromFile(refFileName.c_str());
	} else {
		Log.Message("Encoding reference sequence.");
		std::map<int, RefIdx> binRefMap;

		long size = getSize(Config.GetString("ref"));
		Log.Message("Size of reference genome %ld (%ld)", size, UINT_MAX);
		if (size > UINT_MAX) {
			Log.Error("Reference genome too long! NGM can't handle genomes larger than %ld bytes.", UINT_MAX);
			Fatal();
		}

		uint const binRefSize = ((size / 2) | 1) + 1;
		Log.Message("Allocating %u (%u) bytes for the reference.", binRefSize, FileSize(Config.GetString("ref")));
		binRef = new char[binRefSize];

		gzFile gzfp;
		kseq_t *seq;
		gzfp = gzopen(Config.GetString("ref"), "r");
		seq = kseq_init(gzfp);

		Timer tt;
		tt.ST();
		int l = 0;
		int j = 0;

		char const spacer = 'N';

		//Padding to avoid negative mapping positions
		for (int i = 0; i < 500; ++i) {
			char c = enc4(spacer) << 4;
			c |= enc4(spacer);
			binRef[binRefIndex++] = c;
		}
		while ((l = kseq_read(seq)) >= 0) {
			if(j >= maxRefCount) {
				Log.Error("Currently NextGenMap can't handle more than %d reference sequences.", maxRefCount);
				Fatal();
			}
			binRefMap[j].SeqStart = binRefIndex * 2;
			binRefMap[j].SeqLen = seq->seq.l;
			binRefMap[j].SeqId = j;
			Log.Verbose("Ref %d: %s (%d), Index: %d", j, seq->name.s, seq->seq.l, binRefIndex);
			int nameLength = std::min((size_t) maxRefNameLength, seq->name.l);
			strncpy(binRefMap[j].name, seq->name.s, nameLength);
			binRefMap[j].NameLen = nameLength;
			j += 1;
			char * ref = seq->seq.s;
			for (size_t i = 0; i < seq->seq.l / 2 * 2; i += 2) {
				char c = enc4(ref[i]) << 4;
				c |= enc4(ref[i + 1]);
				binRef[binRefIndex++] = c;
			}
			if (seq->seq.l & 1) {
				char c = enc4(ref[seq->seq.l - 1]) << 4;
				c |= enc4(spacer);
				binRef[binRefIndex++] = c;
			}

			for (int i = 0; i < 500; ++i) {
				//N
				char c = enc4(spacer) << 4;
				c |= enc4(spacer);
				binRef[binRefIndex++] = c;
			}
		}
		refCount = j;
		Log.Message("BinRef length: %d (elapsed %f)", binRefIndex, tt.ET());
		kseq_destroy(seq);
		gzclose(gzfp);

		binRefIndex *= 2;

		if (binRefMap.size() != (size_t) refCount) {
			Log.Error("Error while building ref index.");
			Fatal();
		}
		binRefIdx = new RefIdx[refCount];
		for (int i = 0; i < refCount; ++i) {
			if (binRefMap.find(i) != binRefMap.end()) {
				binRefIdx[i] = binRefMap[i];
//				binRefIdx[i].SeqStart = 0;
				//binRefIdx[i].SeqLen = 2 * binRefIndex - 1;
			} else {
				Log.Error("Error while building ref index.");
				Fatal();
			}
		}
		writeEncRefToFile(refFileName.c_str(), (uint) refCount, binRefSize);
	}

	if (NGM.DualStrand())
	refCount *= 2;

#ifdef VERBOSE
	for (int i = 0; i < refCount; ++i) {
		int len = 0;
		char const * test = GetRefName(i, len);
		Log.Message("%d: Ref: %.*s, Length: %d",i, len, test, (int)GetRefLen(i));
	}
#endif

}

bool _SequenceProvider::DecodeRefSequence(char * const buffer, int n, uint offset, uint bufferLength) {
	uint len = bufferLength - 2;
	if (NGM.DualStrand()) {
		n >>= 1;
	}
//	Log.Message("%u %u %u", offset, bufferLength, binRefIdx[n].SeqLen);
	//if (offset >= binRefIdx[n].SeqLen || offset < 0) {
	if (offset >= GetConcatRefLen() || offset < 0) {
		Log.Verbose("Invalid reference location. Offset: %d", offset);
		return false;
	}
//	int nCount = 0;
//	if (offset < 0) {
//		nCount = abs(offset);
//		len -= nCount;
//		offset = 0;
//	}
	uint end = 0;
	if ((offset + len) > GetConcatRefLen()) {
		end = (offset + len) - GetConcatRefLen();
		len -= end;
	}
//	uint end = std::min((uint)0, binRefIdx[n].SeqLen - (offset + len));
//	if (end < 0) {
//		end = abs(end);
//		len -= end;
//	}
//	int start = binRefIdx[n].SeqStart + ceil(offset / 2.0);
	int start = ceil(offset / 2.0);

	uint codedIndex = 0;
//	for (int i = 0; i < nCount; ++i) {
//		buffer[codedIndex++] = 'x';
//	}
	if (offset & 1) {
		buffer[codedIndex++] = dec4Low(binRef[start - 1]);
	}
	for (uint i = 0; i < ceil(len / 2.0); ++i) {
		buffer[codedIndex++] = dec4High(binRef[start + i]);
		buffer[codedIndex++] = dec4Low(binRef[start + i]);
	}
	if (len & 1) {
		buffer[codedIndex - 1] = 'x';
	}
	for (uint i = 0; i < end; ++i) {
		buffer[codedIndex++] = 'x';
	}

	if (codedIndex > bufferLength) {
		Log.Error("nCount: %d, offset: %d, len: %d (%d), seqlen: %d, end: %d, start: %d, index: %d", 0, offset, bufferLength, (int)ceil(len / 2.0), binRefIdx[n].SeqLen, end, start, codedIndex);
		Log.Error("%.*s", bufferLength, buffer);
		Fatal();
	}

	for (uint i = codedIndex; i < bufferLength; ++i) {
		buffer[i] = '\0';
	}
	return true;
}

char const * _SequenceProvider::GetRefName(int n, int & len) const {
	if (CheckRefNr(n)) {
		if (NGM.DualStrand())
		n >>= 1;
		len = binRefIdx[n].NameLen;

		return binRefIdx[n].name;
	}
	return 0;
}

uint _SequenceProvider::GetConcatRefLen() const {
	return binRefIndex - 1;
}

uint _SequenceProvider::GetRefLen(int n) const {
	if (CheckRefNr(n)) {
		if (NGM.DualStrand())
		n >>= 1;
		return (int) binRefIdx[n].SeqLen;
	} else {
		return 0;
	}
}

uint _SequenceProvider::GetRefStart(int n) const {
	if (CheckRefNr(n)) {
		if (NGM.DualStrand())
		n >>= 1;
		return (int) binRefIdx[n].SeqStart;
	} else {
		return 0;
	}
}

int _SequenceProvider::GetRefCount() const {
	return refCount;
}

bool _SequenceProvider::CheckRefNr(int n) const {
	if (n >= refCount || n < 0) {
		Log.Error("Tried to access invalid reference sequence (%i %x).", n, n);
		throw "Problem";
		//Fatal();
		return false;
	}
	return true;
}

_SequenceProvider::_SequenceProvider() :
		binRef(0), binRefIndex(0), refCount(0), m_EnableBS(false) {
}

_SequenceProvider::~_SequenceProvider() {
	delete[] binRef;
}

void _SequenceProvider::Cleanup() {
	delete pInstance;
}

void _SequenceProvider::PagingUpdate() {

}
