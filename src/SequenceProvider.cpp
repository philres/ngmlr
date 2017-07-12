/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "SequenceProvider.h"

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <stdio.h>
#include <zlib.h>

#include "IConfig.h"
#include "Log.h"
#include "Timing.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#undef module_name
#define module_name "SEQPROV"

_SequenceProvider * _SequenceProvider::pInstance = 0;

_SequenceProvider & _SequenceProvider::Instance() {
	if (pInstance == 0)
		pInstance = new _SequenceProvider();

	return *pInstance;
}

//static int const maxRefCount = 32768;
static int const maxRefCount = 2147483647;
static uint const refEncCookie = 0x74656;

std::string CheckFile(std::string filename, char const * const name) {
	if (!FileExists(filename.c_str())) {
		Log.Error("%s file not found (%s)", name, filename.c_str());
	}
	return filename;
}

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

_SequenceProvider::Chromosome _SequenceProvider::getChrBorders(
		loc start, loc stop) {

	if(start > stop) {
		loc tmp = start;
		start = stop;
		stop = tmp;
	}

	if(start < 1000) {
		start = 1001;
		stop = std::max(1002ll, stop);
	}

//	Log.Message("Detecting chromosome for: %lld %lld", start, stop);

	uloc * upperStart = std::upper_bound(refStartPos, refStartPos + (refCount / 2) + 1, start);
	if ((*upperStart - start) < 1000) {
		upperStart += 1;
//		Log.Message("Start in spacer");
	}
	uloc * upperStop= std::upper_bound(refStartPos, refStartPos + (refCount / 2) + 1, stop);
	// if ((*upperStop - stop) < 1000) don't do anything. Probably still same chromosome

//	Log.Message("Upper: %llu - %llu", *upperStart, *upperStop);

	Chromosome chr;
	if(upperStart == upperStop) {
		chr.start = *(upperStart - 1);
		chr.end = *(upperStart) - 1000;
//		Log.Message("## Chr: %llu - %llu", chr.start, chr.end);
	} else {
//		throw "Could not determine chromosome for interval.";
		chr.start = 0;
		chr.end = 0;
//		Log.Message("## Chr1: %llu - %llu", *(upperStart - 1), *(upperStart) - 1000);
//		Log.Message("## Chr1: %llu - %llu", *(upperStop - 1), *(upperStop) - 1000);
	}

	return chr;
}

_SequenceProvider::Chromosome _SequenceProvider::getChrStart(
		uloc const position) {
	if (position < 1000) {
		Log.Verbose("Can't get starting position of chromosome (%llu).", position);
//		upper += 1;
//		Fatal();
	}
	//Find the next larger chromosome start position in the concatenated reference for the mapping location
	uloc * upper = std::upper_bound(refStartPos,
			refStartPos + (refCount / 2) + 1, position);
	//Check whether the mapping position is in one of the spacer regions between the chromosomes
	if ((*upper - position) < 1000) {
		Log.Verbose("Can't get starting position of chromosome (%llu - %llu = %llu).", *upper, position, *upper - position);
//		Fatal();
		upper += 1;
	}

	Chromosome chr;
	chr.start = *(upper - 1);
	chr.end = *(upper) - 1000;
	return chr;
}

bool _SequenceProvider::convert(SequenceLocation & m_Location) {
	//Convert position back to Chromosome+Position
	SequenceLocation loc = m_Location;
	//Find the next larger chromosome start position in the concatenated reference for the mapping location
	uloc * upper = std::upper_bound(refStartPos,
			refStartPos + (refCount / ((DualStrand) ? 2 : 1)) + 1,
			loc.m_Location);

	//Check whether the mapping position is in one of the spacer regions between the chromosomes
	if ((*upper - loc.m_Location) < 1000) {
		Log.Verbose("Read start position < chromosome start!");
		Log.Verbose("Loc: %u (%d) < %u < %u (%d)", (uloc)*(upper-1), ((upper - 2) - refStartPos) * ((DualStrand) ? 2 : 1), loc.m_Location, (uint)*(upper), ((upper - 1) - refStartPos) * ((DualStrand) ? 2 : 1));

		//Report read/position as unmapped (only happens for --end-to-end)
		return false;
	}

	//Compute actual start position
	loc.m_Location -= *(upper - 1);
	std::ptrdiff_t refId = ((upper - 1) - refStartPos) * ((DualStrand) ? 2 : 1);

	Log.Verbose("Location: %u - Upper: %u - Location: %d, RefId: %d ", m_Location.m_Location, *(upper - 1), loc.m_Location, refId);
	loc.setRefId(refId);
	m_Location = loc;
	return true;
}

int _SequenceProvider::readEncRefFromFile(char const * fileName,
		const uloc maxLen) {
	Log.Message("Reading encoded reference from %s", fileName);
	Timer wtmr;
	wtmr.ST();

	uloc encRefSize = 0;
	uint refCount = 0;
	uint cookie = 0;

	FILE * fp = 0;
	fp = fopen(fileName, "rb");

	size_t read = 0;
	read = fread(&cookie, sizeof(uint), 1, fp);
	read = fread(&refCount, sizeof(uint), 1, fp);
	read = fread(&binRefIndex, sizeof(uloc), 1, fp);
	read = fread(&encRefSize, sizeof(uloc), 1, fp);
	if (cookie != refEncCookie) {
		fclose(fp);
		Log.Error("Invalid encoded reference file found: %s. Please delete it and run NGM again.", fileName);
	}
	if(refCount > maxRefCount) {
		Log.Error("Currently NextGenMap can't handle more than %d reference sequences.", maxRefCount);
	}
	if((encRefSize * 2) > maxLen) {
		Log.Message("Size of reference is %llu Mbp.", encRefSize * 2 / 1000 / 1000);
		Log.Message("With a bin size of 2^%d NextGenMap can only handle a max. reference size of %llu Mbp", Config.getBinSize(), maxLen / 1000 / 1000);
		Log.Message("Please increase --bin-size");
		Log.Error("Max genome size equals 4 GB * 2^bin_size. E.g. with bin size 4 it is 4 GB * 2^4 = 64 GB");
	}
	binRefIdx = new RefIdx[refCount];
	read = fread(binRefIdx, sizeof(RefIdx), refCount, fp);

	binRef = new char[encRefSize];
	read = fread(binRef, sizeof(char), encRefSize, fp);
	fclose(fp);
	Log.Message("Reading %llu Mbp from disk took %.2fs", encRefSize * 2 / 1000 / 1000, wtmr.ET());

	return refCount;
}

void _SequenceProvider::writeEncRefToFile(char const * fileName,
		uint const refCount, uloc const encRefSize) {
	if (!Config.getSkipSave()) {
		Timer wtmr;
		wtmr.ST();
		Log.Message("Writing encoded reference to %s", fileName);
		FILE *fp;

		if (!(fp = fopen(fileName, "wb"))) {
			Log.Warning("WARNING: Could not write encoded reference file to disk: Unable to open output file %s. Please check file permissions.", fileName);
		} else {
			fwrite(&refEncCookie, sizeof(uint), 1, fp);
			fwrite(&refCount, sizeof(uint), 1, fp);
			fwrite(&binRefIndex, sizeof(uloc), 1, fp);
			fwrite(&encRefSize, sizeof(uloc), 1, fp);

			fwrite(binRefIdx, sizeof(RefIdx), refCount, fp);
			fwrite(binRef, sizeof(char), encRefSize, fp);

			fclose(fp);
			Log.Message("Writing to disk took %.2fs", wtmr.ET());
		}
	}
}

uloc getSize(char const * const file) {
	gzFile gzfp;
	kseq_t *seq;
	gzfp = gzopen(Config.getReferenceFile(), "r");
	seq = kseq_init(gzfp);
	int l = 0;
	//1000 -> padding at beginning
	uloc size = 1000;
	while ((l = kseq_read(seq)) >= 0) {
		//1000 -> 1000 x N padding
		int s = (seq->seq.l | 1) + 1;
		size += s + 1000;
	}
	kseq_destroy(seq);
	gzclose(gzfp);
	return size;
}

void _SequenceProvider::Init(bool dualstrand) {
	DualStrand = dualstrand;
	Log.Verbose("Init sequence provider.");

	refFileName = std::string(Config.getReferenceFile()) + std::string("-enc.2.ngm");

	//uint64 used everywhere but in CS rTable, there GetBin division increases range
	const uloc REF_LEN_MAX = UINT_MAX
			* std::max(1.0, pow(2.0, Config.getBinSize()));

	if (FileExists(refFileName.c_str())) {
		//Read
		refCount = readEncRefFromFile(refFileName.c_str(), REF_LEN_MAX);
	} else {
		Log.Message("Encoding reference sequence.");
		std::map<int, RefIdx> binRefMap;

		uloc size = getSize(Config.getReferenceFile());

		Log.Message("Size of reference genome %llu Mbp (max. %llu Mbp)", size / 1000 / 1000, REF_LEN_MAX / 1000 / 1000);

		//Theoretical limit would be uloc max (INT64_MAX), but for speed reasons uint is used in CSTableEntry for bin
		//positions
		if (size > REF_LEN_MAX) {
			Log.Message("With a bin size of 2^%d NextGenMap can only handle a max. reference size of %llu Mbp", Config.getBinSize(), REF_LEN_MAX / 1000 / 1000);
			Log.Message("Please increase --bin-size to handle larger genomes.");
			Log.Error("Max genome size equals 4 GB * 2^bin_size. E.g. with bin size 4 it is 4 GB * 2^4 = 64 GB");
		}

		uloc const binRefSize = ((size / 2) | 1) + 1;
		Log.Verbose("Allocating %llu (%llu) bytes for the reference.", binRefSize, FileSize(Config.getReferenceFile()));
		binRef = new char[binRefSize];

		gzFile gzfp;
		kseq_t *seq;
		gzfp = gzopen(Config.getReferenceFile(), "r");
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
		int skipped = 0;
		while ((l = kseq_read(seq)) >= 0) {
			if(j >= maxRefCount) {
				Log.Error("Currently NextGenMap can't handle more than %d reference sequences.", maxRefCount);
			}
			if(seq->seq.l > minRefSeqLen) {
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
			} else {
				Log.Verbose("Reference sequence %s too short (%d). Skipping.", seq->name.s, seq->seq.l);
				skipped += 1;
			}
		}
		refCount = j;
		Log.Verbose("BinRef length: %ull (elapsed %f)", binRefIndex, tt.ET());
		Log.Message("%d reference sequences were skipped (length < %d).", skipped, minRefSeqLen);
		kseq_destroy(seq);
		gzclose(gzfp);

		binRefIndex = binRefIndex * 2;

		if (binRefMap.size() != (size_t) refCount) {
			Log.Error("Error while building ref index.");
		}
		binRefIdx = new RefIdx[refCount];
		for (int i = 0; i < refCount; ++i) {
			if (binRefMap.find(i) != binRefMap.end()) {
				binRefIdx[i] = binRefMap[i];
//				binRefIdx[i].SeqStart = 0;
				//binRefIdx[i].SeqLen = 2 * binRefIndex - 1;
			} else {
				Log.Error("Error while building ref index.");
			}
		}
		writeEncRefToFile(refFileName.c_str(), (uint) refCount, binRefSize);
	}

	if (DualStrand)
		refCount *= 2;

#ifdef VERBOSE
	for (int i = 0; i < refCount; ++i) {
		int len = 0;
		char const * test = GetRefName(i, len);
		Log.Message("%d: Ref: %.*s, Length: %d",i, len, test, (int)GetRefLen(i));
	}
#endif

	int refCount = SequenceProvider.GetRefCount();
	refStartPos = new uloc[refCount / ((DualStrand) ? 2 : 1) + 1];
	int i = 0;
	int j = 0;
	while (i < refCount) {
		refStartPos[j++] = SequenceProvider.GetRefStart(i);
		i += (DualStrand) ? 2 : 1;
	}
	//Add artificial start position as upper bound for all all reads that map to the last chromosome
	refStartPos[j] = refStartPos[j - 1] + SequenceProvider.GetRefLen(refCount - 1) + 1000;

//	//Test decode
//	char * sequence = new char[10000];
//	int corridor = 12;
//	int decodeLenght = 100;
//
//	DecodeRefSequenceExact(sequence, 2000, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	DecodeRefSequenceExact(sequence, 2001, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	decodeLenght = 101;
//	DecodeRefSequenceExact(sequence, 2000, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	DecodeRefSequenceExact(sequence, 2001, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	DecodeRefSequenceExact(sequence, 1000, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	corridor = 24;
//	DecodeRefSequenceExact(sequence, 1000, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	corridor = 10000;
//	DecodeRefSequenceExact(sequence, 4642603, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	corridor = 10000;
//	DecodeRefSequenceExact(sequence, 1000, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	corridor = 10000;
//	DecodeRefSequenceExact(sequence, 1502, decodeLenght + corridor + 1,
//			corridor);
//	Log.Message("Decoded sequence: %s", sequence);
//
//	exit(0);
}

uloc _SequenceProvider::decode(uloc startPosition, uloc endPosition,
		char * const sequence) {
	Log.Verbose("DEBUG - Start: %llu, End: %llu", startPosition, endPosition);
	uint codedIndex = 0;
	//Position in encoded ref sequence (2 bases per byte)
	uloc start = (startPosition + 1) / 2;//TODO: Check if equiv to ceil(offset / 2.0);
	size_t decodeLength = (endPosition - startPosition + 1) / 2;
	if (startPosition & 1) {
		sequence[codedIndex++] = dec4Low(binRef[start - 1]);
	}
	for (uloc i = 0; i < decodeLength; ++i) {
		sequence[codedIndex++] = dec4High(binRef[start + i]);
		sequence[codedIndex++] = dec4Low(binRef[start + i]);
	}
	return codedIndex;
}

//TODO: remove unnecessary variables
bool _SequenceProvider::DecodeRefSequenceExact(char * const sequence,
		uloc startPosition, uloc sequenceLength, int corridor) {

	if (startPosition >= GetConcatRefLen()) {
		Log.Verbose("Invalid reference location. Offset: %d", startPosition);
		return false;
	}

	memset(sequence, 'x', sequenceLength);

	int halfCorridor = corridor / 2;
	Chromosome chr = getChrStart(startPosition);
	//Log.Message("DEBUG - Chr: %llu - %llu", chr.start, chr.end);

	//First position of the requested sequence (incl corridor)
	uloc decodeStartPosition = startPosition - halfCorridor;

	//Last position of the requested sequence (incl corridor)
	uloc endPosition = startPosition + sequenceLength - halfCorridor;

	uloc decodeEndPosition = endPosition;
	//Requested sequence is longer than the chromosome
	if(endPosition > chr.end) {
		Log.Verbose("Correcting end position");
		uloc diff = endPosition - chr.end;
		//Set end position to last bp of chromosome
		decodeEndPosition -= diff;
	}

	if(halfCorridor > startPosition) {
		Log.Verbose("DecodeStartPosition < 0");
		//DecodeStartPosition < 0
		decodeStartPosition = chr.start;
		uloc diff = halfCorridor - decodeStartPosition + 1000 - (startPosition - chr.start);
		//Log.Message("Start diff: %llu", diff);
		decode(decodeStartPosition, decodeEndPosition, sequence + diff);
	} else if(decodeStartPosition < chr.start) {
		if(decodeEndPosition > chr.start) {
			Log.Verbose("Decoding startposition is in one of the spacer regions");
			//Decoding startposition is in one of the spacer regions
			//Start decoding at chrStartPos, everything before stays N
			uloc diff = chr.start - decodeStartPosition;
			//Log.Message("Start diff: %llu", diff);
			decodeStartPosition += diff;
			decode(decodeStartPosition, decodeEndPosition, sequence + diff);
		} else {
			Log.Verbose("Full interval is in spacer region!");
		}
	} else {
		Log.Verbose("Full decode!");

		SequenceLocation start;
		start.m_Location = decodeStartPosition;

		SequenceLocation end;
		end.m_Location = decodeEndPosition;

		convert(start);
		convert(end);

		int len = 0;

		Log.Verbose("Decoding: %s:%llu-%llu", GetRefName(start.getrefId(), len), start.m_Location, end.m_Location);

		//Decode full sequence
		uloc decodedBp = decode(decodeStartPosition, decodeEndPosition, sequence);
		Log.Verbose("copied %llu to %llu", decodedBp, sequenceLength);
	}

	sequence[sequenceLength - 1] = '\0';

	return true;
}

bool _SequenceProvider::DecodeRefSequence(char * const sequence, int refId,
		uloc position, uloc bufferLength) {
	uloc len = bufferLength - 2;
	if (DualStrand) {
		refId >>= 1;
	}
//	Log.Message("%u %u %u", position, bufferLength, binRefIdx[refId].SeqLen);
	//if (position >= binRefIdx[refId].SeqLen || position < 0) {
	if (position >= GetConcatRefLen()) {
		Log.Verbose("Invalid reference location. Offset: %d", position);
		return false;
	}
//	int nCount = 0;
//	if (position < 0) {
//		nCount = abs(position);
//		len -= nCount;
//		position = 0;
//	}
	uloc end = 0;
	if ((position + len) > GetConcatRefLen()) {
		end = (position + len) - GetConcatRefLen();
		len -= end;
	}
//	uint end = std::min((uint)0, binRefIdx[refId].SeqLen - (position + len));
//	if (end < 0) {
//		end = abs(end);
//		len -= end;
//	}
//	int start = binRefIdx[n].SeqStart + ceil(position / 2.0);
	uloc start = (position + 1) / 2; //TODO: Check if equiv to ceil(position / 2.0);

	uint codedIndex = 0;
//	for (int i = 0; i < nCount; ++i) {
//		sequence[codedIndex++] = 'x';
//	}
	if (position & 1) {
		sequence[codedIndex++] = dec4Low(binRef[start - 1]);
	}
	for (uloc i = 0; i < (len + 1) / 2; ++i) {
		sequence[codedIndex++] = dec4High(binRef[start + i]);
		sequence[codedIndex++] = dec4Low(binRef[start + i]);
	}
	if (len & 1) {
		sequence[codedIndex - 1] = 'x';
	}
	for (uloc i = 0; i < end; ++i) {
		sequence[codedIndex++] = 'x';
	}

	if (codedIndex > bufferLength) {
		Log.Message("nCount: %d, position: %d, len: %d (%d), seqlen: %d, end: %d, start: %d, index: %d", 0, position, bufferLength, (len+1)/2, binRefIdx[refId].SeqLen, end, start, codedIndex);
		Log.Error("%.*s", bufferLength, sequence);
	}

	for (uint i = codedIndex; i < bufferLength; ++i) {
		sequence[i] = '\0';
	}
	return true;
}

char const * _SequenceProvider::GetRefName(int n, int & len) const {
	if (CheckRefNr(n)) {
		if (DualStrand)
			n >>= 1;
		len = binRefIdx[n].NameLen;

		return binRefIdx[n].name;
	}
	return 0;
}

uloc _SequenceProvider::GetConcatRefLen() const {
	return binRefIndex - 1;
}

uloc _SequenceProvider::GetRefLen(int n) const {
	if (CheckRefNr(n)) {
		if (DualStrand)
			n >>= 1;
		return binRefIdx[n].SeqLen;
	} else {
		return 0;
	}
}

uloc _SequenceProvider::GetRefStart(int n) const {
	if (CheckRefNr(n)) {
		if (DualStrand)
			n >>= 1;
		return binRefIdx[n].SeqStart;
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
		return false;
	}
	return true;
}

_SequenceProvider::_SequenceProvider() :
		binRef(0), refStartPos(0), refCount(0) {

	binRefIndex = 0;
	binRefSize = 0;
	binRefIdx = 0;
	DualStrand = true;
}

_SequenceProvider::~_SequenceProvider() {
	delete[] binRef;
	binRef = 0;
	delete[] refStartPos;
	refStartPos = 0;

	if (binRefIdx != 0) {
		delete[] binRefIdx;
		binRefIdx = 0;
	}
}

void _SequenceProvider::Cleanup() {
	delete pInstance;
}

void _SequenceProvider::PagingUpdate() {
	//TODO: remove
}
