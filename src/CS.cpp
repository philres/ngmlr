#include "CS.h"

#include <memory.h>
#include <stdlib.h>
#include <cmath>

#include "NGM.h"
#include "Timing.h"
#include "Debug.h"
#include "AlignmentBuffer.h"

#undef module_name
#define module_name "CS"

//static const int cInvalidLocation = -99999999;
static const SequenceLocation sInvalidLocation(9223372036854775808u, 0, false);

volatile int CS::s_ThreadCount = 0;

bool up = false;
int kCount = 0;

//Default values
int x_SrchTableBitLen = 24;
//Number of bases per batch
int const cBatchSize = 1800000;
//int const cBatchSize =   10000;
uint const cPrefixBaseSkip = 0;

float const cOverflowThresh = 0.1f;
float const cOverflowLimit = 10;

void CS::PrefixIteration(char const * sequence, uloc length,
		PrefixIterationFn func, ulong mutateFrom, ulong mutateTo, void* data,
		uint prefixskip, uloc offset, int prefixBaseCount) {
	prefixBasecount = prefixBaseCount;
	prefixBits = prefixBasecount * 2;
	prefixMask = ((ulong) 1 << prefixBits) - 1;

	CS::PrefixIteration(sequence, length, func, mutateFrom, mutateTo, data,
			prefixskip, offset);

	prefixBasecount = Config.GetInt("kmer", 4, 32);
	prefixBits = prefixBasecount * 2;
	prefixMask = ((ulong) 1 << prefixBits) - 1;
}

// mutate T (0x2) -> C (0x1)
//ulong const mutateFrom = 0x2;
//ulong const mutateTo = 0x1;
//ulong mutateFrom = 0x2;
//ulong mutateTo = 0x1;

void CS::PrefixMutateSearch(ulong prefix, uloc pos, ulong mutateFrom,
		ulong mutateTo, void* data) {
	static int const cMutationLocLimit = Config.GetInt("bs_cutoff");
	ulong const mask = 0x3;

	int mutationLocs = 0;
	for (int i = 0; i < (int) prefixBasecount; ++i) {
		ulong base = mask & (prefix >> (i * 2));
		if (base == mutateFrom)
			++mutationLocs;
	}

	if (mutationLocs <= cMutationLocLimit)
		PrefixMutateSearchEx(prefix, pos, mutateFrom, mutateTo, data);
}

void CS::PrefixMutateSearchEx(ulong prefix, uloc pos, ulong mutateFrom,
		ulong mutateTo, void* data, int mpos) {
	PrefixSearch(prefix, pos, mutateFrom, mutateTo, data);

	ulong const mask = 0x3;
	for (int i = mpos; i < (int) prefixBasecount; ++i) {
		ulong cur = mask & (prefix >> (i * 2));

		if (cur == mutateFrom) {
			ulong p1 = (prefix & ~(mask << (i * 2)));
			ulong p2 = (mutateTo << (i * 2));
			PrefixMutateSearchEx(p1 | p2, pos, mutateFrom, mutateTo, data,
					i + 1);
		}
	}
}

void CS::PrefixSearch(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo,
		void* data) {
	CS * cs = (CS*) data;

	RefEntry const * entries = cs->m_RefProvider->GetRefEntry(prefix,
			cs->m_entry); // Liefert eine liste aller Vorkommen dieses Praefixes in der Referenz
	RefEntry const * cur = entries;

	int const readLength = cs->m_CurrentReadLength;

	if (cur != 0 && cur->refTotal == 0) {
		kCount += 1;
	}

	for (int i = 0; i < cs->m_entryCount; i++) {
		//Get kmer-weight.
		//float weight = cur->weight / 100.0f;
		//cs->weightSum += cur->weight;
		float weight = 1.0f;

		uloc correction =
				(cur->reverse) ?
						(readLength - (pos + CS::prefixBasecount)) : pos;

//#ifdef _DEBUGCSVERBOSE
//		Log.Message("Qry Seq %i - Prefix 0x%x got %i locs (sum %i)", cs->m_CurrentSeq, prefix, cur->refTotal);
//#endif

		int const n = cur->refCount;
//
//		int binNr = pos / 256;
//		int posInBin = pos % 256;
		for (int i = 0; i < n; ++i) {
			uloc loc = cur->getRealLocation(cur->ref[i]);
			//cs->AddLocationStd(GetBin(loc - correction), cur->reverse, weight);
			cs->AddLocationStd(GetBin(loc - (int) (correction * 0.95f)),
					cur->reverse, weight);

//			//Add read pos hashtable
//			uloc key = GetBin(loc - posInBin);
//			if (cur->reverse) {
//				key = key * -1;
//			}
//			if (cs->readBins[binNr].iTable.find(key)
//					== cs->readBins[binNr].iTable.end()) {
//				cs->readBins[binNr].iTable[key] = 0;
//			}
//			cs->readBins[binNr].iTable[key] += 1;
//			cs->readBins[binNr].n += 1;
//			cs->readBins[binNr].max = std::max(cs->readBins[binNr].max,
//					cs->readBins[binNr].iTable[key]);
		}

		cur++;
	}
}

void CS::AddLocationStd(uloc const m_Location, bool const reverse,
		double const freq) {
	bool newEntry = false;

	CSTableEntry* entry = rTable + Hash(m_Location);
	CSTableEntry* const maxEntry = rTable + c_SrchTableLen;

	//uint hpo = Hash( m_Location );
	while ((newEntry = (entry->state & 0x7FFFFFFF) == currentState)
			&& !(entry->m_Location == m_Location)) {
		++entry;
		if (entry >= maxEntry)
			entry = rTable;
		if (--hpoc == 0) {
			throw 1;
		}
	}

	float score = freq;
	if (!newEntry) {
		entry->m_Location = m_Location;
		entry->state = currentState & 0x7FFFFFFF;
		if (reverse) {
			entry->fScore = 0.0f;
			entry->rScore = score;
		} else {
			entry->fScore = score;
			entry->rScore = 0.0f;
		}
	} else {
		//Add kmer-weight to position
		if (reverse) {
			score = (entry->rScore += freq);
		} else {
			score = (entry->fScore += freq);
		}
	}

	//Compute max k-mer weight for this read
	if (score > maxHitNumber) {
		maxHitNumber = score;
		//currentThresh = round(maxHitNumber * m_CsSensitivity);
		currentThresh = maxHitNumber * m_CsSensitivity;
	}

	//If kmer-weight larger than threshold -> add to rList.
	if (!(entry->state & 0x80000000) && score >= currentThresh) {
		entry->state |= 0x80000000;
		rList[rListLength++] = entry - rTable;
	}

}

void CS::AddLocationFallback(SequenceLocation const & loc, double const freq) {
	throw "Not implemented (AddLocationFallback)";
}

#include <stdio.h>
#include <sstream>

void CS::debugCS(MappedRead * read, int& n, float& mi_Threshhold) {

	std::stringstream ss;
	ss << Config.GetString("output") << "_kmer-profile/" << read->ReadId
			<< ".csv";
	Log.Message("Opening %s", ss.str().c_str());
	FILE * kmerCount = fopen(ss.str().c_str(), "w");

	int count = 0;
	int accepted = 0;
	float max = 0.0f;
	for (int i = 0; i < c_SrchTableLen; ++i) {

		if ((rTable[i].state & 0x7FFFFFFF) == currentState) {

			max = std::max(std::max(max, rTable[i].fScore), rTable[i].rScore);
			//Correct position
			SequenceLocation loc;
			loc.m_Location = ResolveBin(rTable[i].m_Location);
//			Log.Message("Could not convert location %llu", ResolveBin(rTable[i].m_Location));
//			if (!SequenceProvider.convert(loc)) {
//				Log.Message("ERROR");
//				loc.m_Location = ResolveBin(rTable[i].m_Location);
//				loc.setRefId(0);
//
//			}

			int refNameLength = 0;
			if (rTable[i].fScore > 0.0f) {
				if (rTable[i].fScore >= mi_Threshhold) {
					Log.Debug(8192, "READ_%d\tCS_RESULTS\tInternal location: %llu (+), Location: %llu (Ref: %s), Score: %f (ACCEPT)", read->ReadId, rTable[i].m_Location, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), rTable[i].fScore);
//					Log.Message("READ_%d\tCS_RESULTS\tInternal location: %llu (+), Location: %llu (Ref: %s), Score: %f (ACCEPT)", read->ReadId, rTable[i].m_Location, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), rTable[i].fScore);
					accepted += 1;
					fprintf(kmerCount, "%s;%llu;%f;1;%d\n", read->name, ResolveBin(rTable[i].m_Location), rTable[i].fScore, read->length);
				} else {
					Log.Debug(4096, "READ_%d\tCS_DETAILS\tInternal location: %llu (+), Location: %llu (Ref: %s), Score: %f (REJECT)", read->ReadId, rTable[i].m_Location, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), rTable[i].fScore);
					fprintf(kmerCount, "%s;%llu;%f;0;%d\n", read->name, ResolveBin(rTable[i].m_Location), rTable[i].fScore, read->length);
				}
				count += 1;
			}
			if (rTable[i].rScore > 0.0f) {
				if (rTable[i].rScore >= mi_Threshhold) {
					Log.Debug(8192, "READ_%d\tCS_RESULTS\tInternal location: %llu (-), Location: %llu (Ref: %s), Score: %f (ACCEPT)", read->ReadId, rTable[i].m_Location, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), rTable[i].rScore * -1.0f);
//					Log.Message("READ_%d\tCS_RESULTS\tInternal location: %llu (-), Location: %llu (Ref: %s), Score: %f (ACCEPT)", read->ReadId, rTable[i].m_Location, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), rTable[i].rScore);
					accepted += 1;
					fprintf(kmerCount, "%s;%llu;%f;1;%d\n", read->name, ResolveBin(rTable[i].m_Location), rTable[i].rScore * -1.0f, read->length);
				} else {
					Log.Debug(4096, "READ_%d\tCS_DETAILS\tInternal location: %llu (-), Location: %llu (Ref: %s), Score: %f (REJECT)", read->ReadId, rTable[i].m_Location, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), rTable[i].rScore * -1.0f);
					fprintf(kmerCount, "%s;%llu;%f;0;%d\n", read->name, ResolveBin(rTable[i].m_Location), rTable[i].rScore * -1.0f, read->length);
				}
				count += 1;
			}
		}
	}

	Log.Debug(8, "%s\tCS\tSummary: %d CMRs (%d without filter), %d accepted, theta: %f, max: %f", read->name, n, count, accepted, mi_Threshhold, max);
	fclose(kmerCount);
}

int CS::CollectResultsStd(MappedRead * read) {

	static float const mink = Config.GetFloat("kmer_min");

	int max = (read->length - CS::prefixBasecount + 1) * 0.9;

	if (kCount > max)
		read->mappingQlty = 0;

	read->s = maxHitNumber;

	float mi_Threshhold = std::max(mink, currentThresh);

	int n = rListLength;

//#ifdef DEBUGLOG
	//if(Config.GetInt(LOG_LVL) > 0) {
//if (strcmp(read->name, "SRR1284073.163322") == 0) {
	//debugCS(read, n, mi_Threshhold);
//}
	//}
//#endif

	int index = 0;
	if (2 * n > tmpSize) {
		tmpSize = 2 * n * 1.5f;
		Log.Debug(LOG_CS_DETAILS, "Increasing size of location buffer to %d", tmpSize);
		delete[] tmp;
		tmp = new LocationScore[tmpSize];
	}

	for (int i = 0; i < n; ++i) {
		CSTableEntry temp = rTable[rList[i]];

		if (temp.fScore >= mi_Threshhold) {
			LocationScore * toInsert = &tmp[index++];
			toInsert->Score.f = temp.fScore;
//			toInsert->Location.m_Location = ResolveBin(temp.m_Location);
			toInsert->Location.m_Location = temp.m_Location;
			toInsert->Location.setReverse(false);
		}
		if (temp.rScore >= mi_Threshhold) {
			LocationScore * toInsert = &tmp[index++];
			toInsert->Score.f = temp.rScore * -1.0f;
//			toInsert->Location.m_Location = ResolveBin(temp.m_Location);
			toInsert->Location.m_Location = temp.m_Location;
			toInsert->Location.setReverse(true);
		}
	}
	static int const maxScores = Config.GetInt("max_cmrs");
	if (index < maxScores)
		read->AllocScores(tmp, index);

	return index;
}

int CS::CollectResultsFallback(MappedRead * read) {
	throw "Not implemented (CollectResultsFallback)";
	return 0;
}

bool sortLocationScore2(LocationScore a, LocationScore b) {
	return a.Location.m_Location < b.Location.m_Location;
}

bool sortLocationScoreByScore(LocationScore a, LocationScore b) {
	return fabs(a.Score.f) > fabs(b.Score.f);
}

template<typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

Align CS::computeAlignment(MappedRead* read, int const scoreId,
		int const corridor) {

	Align align;
	LocationScore & score = read->Scores[scoreId];

	Log.Message("Computing alignment (%d) for position: %llu", scoreId, score.Location.m_Location);
	char * refBuffer = new char[read->length + corridor + 1];

	//decode reference sequence
	if (!SequenceProvider.DecodeRefSequence(refBuffer, 0,
			score.Location.m_Location - (corridor >> 1), read->length + corridor)) {
		//Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
		Log.Warning("Could not decode reference for alignment (read: %s)", read->name);
		memset(refBuffer, 'N', read->length * 1.2f);
	}
	//initialize arrays for CIGAR and MD string
	align.pBuffer1 = new char[read->length * 4];
	align.pBuffer2 = new char[read->length * 4];
	*(int*) align.pBuffer1 = 0x212121;
	*(int*) align.pBuffer2 = 0x212121;

	int const mode = 0;
	if (score.Location.isReverse()) {
		oclAligner->SingleAlign(mode, corridor, (char const * const ) refBuffer,
				(char const * const ) read->RevSeq, align, 0);
	} else {
		oclAligner->SingleAlign(mode, corridor, (char const * const ) refBuffer,
				(char const * const ) read->Seq, align, 0);
	}

	score.Location.m_Location += align.PositionOffset - (corridor >> 1);

	delete[] refBuffer;

	return align;
}

struct interval {
	int start;
	int end;
};

void CS::BuildReducedRef(ulong prefix, uloc pos, ulong mutateFrom,
		ulong mutateTo, void* data) {

	CS * cs = (CS*) data;

	if (cs->reducedRefTab[prefix].n != 10000) {
		cs->reducedRefTab[prefix].positions[cs->reducedRefTab[prefix].n++] = pos
				+ cs->currentReducedRefOffset;
		//reducedRefTab[0].positions[0] = pos;
	}
	//Log.Message("%llu at %llu", prefix, uloc);

}

static const unsigned char ReverseTable16rt[] = { 0x00, 0x04, 0x08, 0x0C, 0x01,
		0x05, 0x09, 0x0D, 0x02, 0x06, 0x0A, 0x0E, 0x03, 0x07, 0x0B, 0x0F };

//Works only for 4 byte
inline ulong revCompRt(ulong prefix) {
	static const int shift = 32 - CS::prefixBits;

	//Compute complement
	ulong compPrefix = (prefix ^ 0xAAAAAAAA) & CS::prefixMask;
	//Reverse
	compPrefix = compPrefix << shift;
	ulong compRevPrefix = (ReverseTable16rt[compPrefix & 0x0f] << 28)
			| (ReverseTable16rt[(compPrefix >> 4) & 0x0f] << 24)
			| (ReverseTable16rt[(compPrefix >> 8) & 0x0f] << 20)
			| (ReverseTable16rt[(compPrefix >> 12) & 0x0f] << 16)
			| (ReverseTable16rt[(compPrefix >> 16) & 0x0f] << 12)
			| (ReverseTable16rt[(compPrefix >> 20) & 0x0f] << 8)
			| (ReverseTable16rt[(compPrefix >> 24) & 0x0f] << 4)
			| (ReverseTable16rt[(compPrefix >> 28) & 0x0f]);

//	ulong compRevPrefix = (ReverseTable256[compPrefix & 0xff] << 24)
//			| (ReverseTable256[(compPrefix >> 8) & 0xff] << 16)
//			| (ReverseTable256[(compPrefix >> 16) & 0xff] << 8)
//			| (ReverseTable256[(compPrefix >> 24) & 0xff]);

	return compRevPrefix;
}

void CS::AddLocationRead(ulong prefix, uloc pos, ulong mutateFrom,
		ulong mutateTo, void* data) {
	CS * cs = (CS*) data;

	int binNr = pos / 256;
	int posInBin = pos % 256;

	posInBin = posInBin * 0.95f;

	int kmerCount = cs->reducedRefTab[prefix].n;

	for (int i = 0; i < kmerCount; ++i) {
		uloc loc = cs->reducedRefTab[prefix].positions[i];
		uloc key = (loc - posInBin) / cs->reducedBinSize;
		if (cs->readBins[binNr].iTable.find(key)
				== cs->readBins[binNr].iTable.end()) {
			cs->readBins[binNr].iTable[key].count = 0;
			cs->readBins[binNr].iTable[key].sum = 0;
		}
		cs->readBins[binNr].iTable[key].count += 1;
		cs->readBins[binNr].iTable[key].sum += (loc - posInBin);
		cs->readBins[binNr].n += 1;
		cs->readBins[binNr].max = std::max(cs->readBins[binNr].max,
				cs->readBins[binNr].iTable[key].count);
	}

	ulong revPrefix = revCompRt(prefix);
	char * seq = cs->m_CurrentBatch[cs->m_CurrentSeq]->Seq;

//	Log.Message("Seq: %.6s, Fwd: %llu, Rev: %llu", seq + pos, prefix, revPrefix);
	kmerCount = cs->reducedRefTab[revPrefix].n;

	for (int i = 0; i < kmerCount; ++i) {
		uloc loc = cs->reducedRefTab[revPrefix].positions[i];
//		uloc key = (loc - posInBin) / cs->reducedBinSize;

//		Log.Message("%d-%d - correction : %d", binNr, posInBin, (256 - (posInBin + CS::prefixBasecount)));
		uloc key = (loc - (256 - (posInBin + CS::prefixBasecount)))
				/ cs->reducedBinSize;
		key = key * -1;
//		Log.Message("%llu -> %llu", loc, key);
		if (cs->readBins[binNr].iTable.find(key)
				== cs->readBins[binNr].iTable.end()) {
			cs->readBins[binNr].iTable[key].count = 0;
			cs->readBins[binNr].iTable[key].sum = 0;
		}
		cs->readBins[binNr].iTable[key].count += 1;
		cs->readBins[binNr].iTable[key].sum += (loc - (256 - (posInBin + CS::prefixBasecount)));
		cs->readBins[binNr].n += 1;
		cs->readBins[binNr].max = std::max(cs->readBins[binNr].max,
				cs->readBins[binNr].iTable[key].count);
	}
}

void CS::SendToBuffer(MappedRead * read, ScoreBuffer * sw,
		AlignmentBuffer * out) {

	float maxScore = 0.0f;

	int sumScoreIndex = 0;

	if (read->numScores() > 0) {
		LocationScore * tmpScores = new LocationScore[read->numScores()];
		memcpy(tmpScores, read->Scores,
				sizeof(LocationScore) * read->numScores());

		//TODO: sort by location
		std::sort(tmpScores, tmpScores + read->numScores(), sortLocationScore2);

		for (int i = 0; i < read->numScores(); ++i) {

			LocationScore score = tmpScores[i];

			SequenceLocation loc = score.Location;
			//SequenceProvider.convert(loc);

			int refNameLength = 0;
			Log.Message("READ_%d\tSCORES_RESULTS\tCMR_%d\t%f\t%llu\t%s", read->ReadId, i, score.Score.f, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength));
		}

		//Merge peaks
		for (int i = 0; i < read->numScores(); ++i) {

			uloc position = tmpScores[i].Location.m_Location;
			LocationScore score = tmpScores[i];
			float scoresum = tmpScores[i].Score.f;

			int j = i + 1;
			while (j < read->numScores()
					&& abs(
							tmpScores[i].Location.m_Location
									- tmpScores[j].Location.m_Location) < 3
					&& sgn(tmpScores[i].Score.f) == sgn(tmpScores[j].Score.f)) {
				scoresum += tmpScores[j].Score.f;
				j += 1;
				i += 1;
			}

			SequenceLocation loc = score.Location;
			loc.m_Location = ResolveBin(loc.m_Location);
//			SequenceProvider.convert(loc);

			read->Scores[sumScoreIndex].Location = loc;
			read->Scores[sumScoreIndex].Score.f = scoresum;

			maxScore = std::max(maxScore, (float) fabs(scoresum));

			Log.Message("%llu: %f", read->Scores[sumScoreIndex].Location.m_Location, read->Scores[sumScoreIndex].Score.f);

			sumScoreIndex += 1;

		}
//		printf("%s\t%d\t%f\t%d\t%f\n", read->name, count, maxScore,
//				read->length, maxScore / read->length);
	}

	Log.Message("Waiting");
	getchar();

	uloc reducedRefStart = 0;
	uloc reducedRefStop = 0;
	Log.Message("Sum score index: %d", sumScoreIndex);
	//CMR serach on reduced reference with smaller k-mers
	for (int i = 0; i < sumScoreIndex; ++i) {
		Log.Message("%f > %f", read->Scores[i].Score.f, maxScore * 0.99f);
		if (abs(read->Scores[i].Score.f) > maxScore * 0.99f) {
			Log.Message("Candidate %llu: %f", read->Scores[i].Location.m_Location, read->Scores[i].Score.f);
//			reducedRefStart = read->Scores[i].Location.m_Location - read->length / 2;
//			reducedRefStop = reducedRefStart + 2 * read->length;
			reducedRefStart = read->Scores[i].Location.m_Location - 50000;
			reducedRefStop = reducedRefStart + 2 * 50000;
			Log.Message("Reduced reference from %llu to %llu", reducedRefStart, reducedRefStop);
		}
	}

	reducedBinSize = 16;
	int maxKmerPos = 10000;
	double kmerSize = 6;
	int reducedRefTabLen = pow(4.0, kmerSize);
	reducedRefTab = new RefTabEntry[reducedRefTabLen + 1];
	currentReducedRefOffset = reducedRefStart;

	//Init reduced ref tab
	for (int k = 0; k < reducedRefTabLen; ++k) {
		reducedRefTab[k].positions = new int[maxKmerPos];
		reducedRefTab[k].n = 0;
	}

	Timer t;
	t.ST();
	Log.Message("Start timing");
	int reducedRefTabSeqLen = reducedRefStop - reducedRefStart + 1;
	char * reducedRefTabSeq = new char[reducedRefTabSeqLen];

	SequenceProvider.DecodeRefSequence(reducedRefTabSeq,0, reducedRefStart, reducedRefTabSeqLen);
	//Log.Message("%s", reducedRefTabSeq);

	//Building hash table
	CS::prefixBasecount = kmerSize;
	CS::prefixBits = prefixBasecount * 2;
	CS::prefixMask = ((ulong) 1 << prefixBits) - 1;

	Timer t1;
	t1.ST();
	PrefixIteration(reducedRefTabSeq, reducedRefTabSeqLen, BuildReducedRef, 0,
			0, this, 0, 0);
	Log.Message("Build reduced ref table took: %f", t1.ET());

//	long sum = 0;
//	for (int i = 0; i < reducedRefTabLen; ++i) {
//		sum += reducedRefTab[i].n;
//	}
//	Log.Message("Mean: %f", sum / reducedRefTabLen * 1.0f);

	Timer t2;
	t2.ST();
	PrefixIteration(read->Seq, read->length, AddLocationRead, 0, 0, this, 0, 0);
	Log.Message("Adding locations took %f seconds", t2.ET());

	CS::prefixBasecount = 13;
	CS::prefixBits = prefixBasecount * 2;
	CS::prefixMask = ((ulong) 1 << prefixBits) - 1;

	Log.Message("Hash loop up finished");

	//Finding positions
	Timer t3;
	t3.ST();
	int binCount = read->length / 256 + 1;
	LocationScore * scores = new LocationScore[100];
	int scoreIndex = 0;

	uloc lastBinPos = 0;

	for (int j = 0; j < binCount; ++j) {
		scoreIndex = 0;
		ReadBin & bin = readBins[j];
		//Log.Message("Bin %d: %d, %d (elements in hash: %d)", j, bin.n, bin.max, bin.iTable.size());
		if (bin.n == -1) {
			Log.Message("Overflow");
		} else {
			if(bin.n > 0) {
				int k = 0;

				//Find highest scoring positions
				typedef std::map<long, AvgPos>::iterator it_type;
				for(it_type iterator = bin.iTable.begin(); iterator != bin.iTable.end(); iterator++) {
					int count = iterator->second.count;
					uloc avgPos = iterator->second.sum / count;
					//If the position has at least half as many k-mer hits as the best position
					if(count > bin.max * 0.5f) {

						scores[scoreIndex].Score.i = count;
						scores[scoreIndex].Location.m_Location = iterator->first * reducedBinSize;

						//float alnScore = scoreReadPart(read->Seq, read->length, j, iterator->first * reducedBinSize);
						float alnScore = 0;

						SequenceLocation loc;
						loc.m_Location = abs(avgPos);
						if(avgPos < 0) {
							loc.setReverse(true);
						}
						//loc.m_Location = abs(iterator->first * reducedBinSize);
						//SequenceProvider.convert(loc);
						//int refNameLength = 0;

//						if(count == bin.max) {
//							Log.Message("BIN %d - %d: from %d to %d (blocknr: %d) x %d (max: %d) - score: %f - %s:%llu", j, k++, iterator->first * reducedBinSize, iterator->first * reducedBinSize + 256, iterator->first * reducedBinSize / 256, count, bin.max, alnScore, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), loc.m_Location);
//							Log.Message("Start: %llu", avgPos);
//						} else {
//							Log.Message("BIN %d - %d: \tfrom %d to %d (blocknr: %d) x %d (max: %d) - score: %f - %s:%llu", j, k++, iterator->first * reducedBinSize, iterator->first * reducedBinSize + 256, iterator->first * reducedBinSize / 256, count, bin.max, alnScore, SequenceProvider.GetRefName(loc.getrefId(), refNameLength), loc.m_Location);
//							Log.Message("Start: %llu", avgPos);
//						}

						scoreIndex += 1;
					} else {
						//Log.Message("BIN %d - %d: \t\t\t%llu x %d (max: %d)", j, k++, iterator->first, count, bin.max);
					}
				}

				k = 0;
				int bestIndex = -1;
				if(scoreIndex == 1) {
					bestIndex = 0;
					k += 1;
				} else {
					float bestScore = -1;
					for(int l = 0; l < scoreIndex; l++) {
						float alnScore = scoreReadPart(read->Seq, read->length, j, scores[l].Location.m_Location);
						Log.Message("\tSCORE: %llu %d -> %f", scores[l].Location.m_Location, scores[l].Score.i, alnScore);
						k += 1;
						if(alnScore > bestScore) {
							bestIndex = l;
							bestScore = alnScore;
						}
					}

				}
				if(k == 0) {
					Log.Message("BIN %d - iTable was empty", j);
				}

				Log.Message("BIN %d - %d: %d, %d (blocknr: %d) x %d (max: %d)", j, k++, scores[bestIndex].Location.m_Location, scores[bestIndex].Location.isReverse(), scores[bestIndex].Location.m_Location / 256, scores[bestIndex].Score.i, bin.max);

				if(j == 0) {
					//First BIN
					Log.Message("Starting at %llu", scores[bestIndex].Location.m_Location);
				} else {
					uloc offset = scores[bestIndex].Location.m_Location - lastBinPos;
					Log.Message("offset %d", offset);
					if(offset > 300 || offset < 200) {
						SequenceLocation conv;
						conv.m_Location = lastBinPos;
						SequenceProvider.convert(conv);
						int length = 0;
						Log.Message("Offset ist off at position %s-%llu", SequenceProvider.GetRefName(conv.getrefId(), length), conv.m_Location);
						getchar();
					}

				}

				lastBinPos = scores[bestIndex].Location.m_Location;
			} else {
				Log.Message("Bin.n == 0");
			}
		}
		Log.Message("----------------------------------------------------");
	}
	delete[] scores;
	scores = 0;
	Log.Message("Assigning bins to reference took %f", t3.ET());

	Log.Message("Took %f secdons", t.ET());

	//Delete reduced ref tab
	for (int k = 0; k < reducedRefTabLen; ++k) {
		delete[] reducedRefTab[k].positions;
		reducedRefTab[k].positions = 0;
	}
	delete[] reducedRefTab;
	reducedRefTab = 0;

	delete[] reducedRefTabSeq;
	reducedRefTabSeq = 0;

	Log.Message("Waiting");
	getchar();

	if (sumScoreIndex < 1) {
		Log.Message("Read %s has no candidate", read->name);
		alignmentBuffer->WriteRead(read, false);
	} else if (sumScoreIndex < 100) {
		Log.Message("Read %s has %d candidates", read->name, sumScoreIndex);
		std::sort(read->Scores, read->Scores + sumScoreIndex,
				sortLocationScoreByScore);

		read->Calculated = sumScoreIndex;
		read->Alignments = new Align[read->Calculated];
		read->numTopScores = 1;

		read->computeReverseSeq();

		int corridor = (int) (read->length * 0.2f + 0.5f) / 2 * 2;

		bool finished = false;
		for (int i = 0; i < sumScoreIndex && !finished; ++i) {
			int refNameLength = 0;
//			printf("%s\t%d\t%s\t%llu\t%f\n", read->name, read->length,
//					SequenceProvider.GetRefName(read->Scores[i].Location.getrefId(), refNameLength), read->Scores[i].Location.m_Location, read->Scores[i].Score.f);

			float kmerPer100Bp = read->Scores[i].Score.f * 100.0f / read->length;
			Log.Message("Found candidate (%d) with k-mer score %f (%f) on position %llu", i, read->Scores[i].Score.f, kmerPer100Bp, read->Scores[i].Location.m_Location);

			Align align = computeAlignment(read, i, corridor);
			read->Alignments[i] = align;

//			SequenceProvider.convert(read->Scores[i].Location);

			int clipped = align.QStart + align.QEnd;
			float covered = (read->length - clipped) * 100.0f / read->length;

			float scorePer100Bp = read->Alignments[i].Score * 100.0f / (read->length - clipped);
			Log.Message("Alignment with score %f (%f) to position %llu (%d) covered %f%% of the read (%d: %d - %d)", align.Score, scorePer100Bp, read->Scores[i].Location.m_Location, !read->Scores[i].Location.isReverse(), covered, read->length - clipped, align.QStart, read->length - align.QEnd);
			//finished = true;
			if(covered > 90.0) {
				finished = true;
			}
		}

		alignmentBuffer->WriteRead(read, true);
//		NGM.AddMappedRead(read->ReadId);
	} else {

		Log.Message("Read %s has to many candidates (%d)", read->name, sumScoreIndex);
		alignmentBuffer->WriteRead(read, false);
		//NGM.AddUnmappedRead(read, 0);
		//delete read;
		//read = 0;

	}

	return;

//	if (read == 0)
//		return;

//	int count = read->numScores();
//
//	if (count == 0) {
//		read->Calculated = 0;
//		out->addRead(read, -1);
//	} else {
//		read->Calculated = 0;
//		sw->addRead(read, count);
//		++m_WrittenReads;
//	}
}

void CS::Cleanup() {
	NGM.ReleaseWriter();
}

float CS::scoreReadPart(char const * const readSeq, int const qryLen,
		int const bin, long const refBin) {

	char * readPart = new char[256 + 1];
	//Log.Message("Copy %d bytes starting from +%d", 256, bin * 256);
	strncpy(readPart, readSeq + bin * 256, 256);
	char * refPart = new char[256 + 16 + 1];
	SequenceProvider.DecodeRefSequence(refPart, 0, refBin - 8, 256 + 16);

//	SequenceLocation loc;
//	loc.m_Location = refBin - 8;
//	SequenceProvider.convert(loc);

	//Log.Message("Ref from: %d to %d", loc.m_Location, loc.m_Location + 256 + 16);

	//Log.Message("ReadBin: %d, RefBin: %d", bin, refBin);
	//Log.Message("Ref:  %s\nRead: %s", refPart, readPart);

	float score = 0;
	oclAligner->SingleScore(0, reducedBinSize, refPart, readPart, score, 0);

	delete[] readPart;
	delete[] refPart;

	return score;
}

int CS::RunBatch(ScoreBuffer * sw, AlignmentBuffer * out) {
	PrefixIterationFn pFunc = &CS::PrefixSearch;
	if (m_EnableBS)
		pFunc = &CS::PrefixMutateSearch;

	int nScoresSum = 0;

	for (size_t i = 0; i < m_CurrentBatch.size(); ++i) {
		Log.Message("Processing read %s", m_CurrentBatch[i]->name);
		int binCount = m_CurrentBatch[i]->length / 256 + 1;
		//readBins = new std::vector<int>[binCount + 1];
		readBins = new ReadBin[binCount + 1];
		for (int j = 0; j < binCount; ++j) {
			ReadBin * bin = readBins + j;
			//bin->slots = new uloc[100000];
			bin->max = 0;
			bin->n = 0;
		}
		Log.Message("BinCount: %d (read length: %d)", binCount, m_CurrentBatch[i]->length);

//		m_CurrentSeq = m_CurrentBatch[i]->ReadId;
		m_CurrentSeq = i;

		//currentState = 2 * i;
		currentState++;
		rListLength = 0;
		maxHitNumber = 0.0f;
		currentThresh = 0.0f;
		//weightSum = 0.0f;
		kCount = 0;

		ulong mutateFrom;
		ulong mutateTo;
		static bool const isPaired = Config.GetInt("paired") > 0;
		if (isPaired && (m_CurrentBatch[i]->ReadId & 1)) {
			//Second mate
			mutateFrom = 0x0;
			mutateTo = 0x3;
		} else {
			//First mate
			mutateFrom = 0x2;
			mutateTo = 0x1;
		}

		m_CurrentReadLength = m_CurrentBatch[i]->length;

		++m_ProcessedReads;
		bool fallback = m_Fallback;

		char const * const qrySeq = m_CurrentBatch[i]->Seq;
		uloc qryLen = m_CurrentBatch[i]->length;

		SetSearchTableBitLen(24);

		if (!fallback) {
			try {
				hpoc = c_SrchTableLen * 0.333f;
				PrefixIteration(qrySeq, qryLen, pFunc, mutateFrom, mutateTo,
						this, m_PrefixBaseSkip);
				nScoresSum += CollectResultsStd(m_CurrentBatch[i]);
			} catch (int overflow) {
				++m_Overflows;

				Log.Message("Overflow in CS!");
				// Initiate fallback
				fallback = true;
			}
		}
		if (fallback) {
			int c_SrchTableBitLenBackup = c_SrchTableBitLen;
			int x = 2;
			while (fallback && (c_SrchTableBitLenBackup + x) <= 24) {
				fallback = false;
				try {

					SetSearchTableBitLen(c_SrchTableBitLenBackup + x);

					rListLength = 0;
					//currentState = 2 * i + 1;
					currentState++;
					maxHitNumber = 0.0f;
					currentThresh = 0.0f;

					hpoc = c_SrchTableLen * 0.777f;
					PrefixIteration(qrySeq, qryLen, pFunc, mutateFrom, mutateTo,
							this, m_PrefixBaseSkip);
					nScoresSum += CollectResultsStd(m_CurrentBatch[i]);

				} catch (int overflow) {
					Log.Message("Overflow in fallback!");
					fallback = true;
					x += 1;
				}
			}
			SetSearchTableBitLen(c_SrchTableBitLenBackup);
		}

//		for (int j = 0; j < binCount; ++j) {
//
//			ReadBin & bin = readBins[j];
//			if(bin.n == -1) {
//				Log.Message("Overflow");
//			} else {
//				if(bin.n > 0) {
//					std::sort(bin.slots, bin.slots + bin.n);
//					int count = 1;
//					int last = -1;
//					for(int k = 0; k < bin.n; ++k) {
//						if(last == bin.slots[k]) {
//							count += 1;
//						} else {
//							if(count > 5) {
//								Log.Message("%d: %llu x %d", j, last, count);
//							}
//							last = bin.slots[k];
//							count = 1;
//						}
//					}
//					Log.Message("%d: %llu x %d", j, last, count);
//
//					delete readBins[j].slots;
//					readBins[j].slots = 0;
//				}
//			}
//		}

		SendToBuffer(m_CurrentBatch[i], sw, out);

		delete[] readBins;
		readBins = 0;
		Log.Message("Waiting");
		getchar();

	}

	return nScoresSum;
}

void CS::DoRun() {

	int gpu = m_TID;
	if (Config.Exists("gpu")) {
		int threadcount = 1;
		static const int cMaxAligner = 32;
		threadcount = Config.GetInt("gpu", 1, cMaxAligner);
		int * gpus = new int[threadcount];
		if (Config.Exists("gpu")) {
			Config.GetIntArray("gpu", gpus, threadcount);
		} else {
			gpus[0] = 0;
		}
		gpu = gpus[m_TID % threadcount];
	}

	NGM.AquireOutputLock();
	oclAligner = NGM.CreateAlignment(gpu | (std::min(Config.GetInt("format", 0, 2), 1) << 8));
	alignmentBuffer = new AlignmentBuffer(
	Config.Exists("output") ? Config.GetString("output") : 0, oclAligner);
	ScoreBuffer * scoreBuffer = new ScoreBuffer(oclAligner, alignmentBuffer);
	NGM.ReleaseOutputLock();

	int x_SrchTableLen = (int) pow(2, x_SrchTableBitLen);

	rTable = new CSTableEntry[x_SrchTableLen];
	Log.Debug(LOG_CS_DETAILS, "Sizeof CSTableEntry %d (%d)", sizeof(CSTableEntry), sizeof(SequenceLocation));
	Log.Debug(LOG_CS_DETAILS, "rTable: %d (%d x (%d + %d))", (sizeof(CSTableEntry) + sizeof(int)) * x_SrchTableLen, x_SrchTableLen, sizeof(CSTableEntry), sizeof(int));
	rList = new int[x_SrchTableLen];

	for (int i = 0; i < x_SrchTableLen; ++i) {
		rTable[i].m_Location = sInvalidLocation.m_Location;
		rTable[i].state = -1;
		rList[i] = -1;
	}

	m_CsSensitivity = 0.50f;
	if (Config.Exists("sensitivity")) {
		m_CsSensitivity = Config.GetFloat("sensitivity", 0, 1);
	} else {
		Log.Warning("Sensitivity parameter neither set nor estimated. Falling back to default.");
	}

	m_RefProvider = 0;

	Timer tmr;
	m_RefProvider = NGM.GetRefProvider(m_TID);
	AllocRefEntryChain();
	while (NGM.ThreadActive(m_TID, GetStage()) && ((m_CurrentBatch = NGM.GetNextReadBatch(m_BatchSize)), (m_CurrentBatch.size() > 0))) {
		Log.Debug(LOG_CS_DETAILS, "CS Thread %i got batch (len %i)", m_TID, m_CurrentBatch.size());
		tmr.ST();

		int nCRMsSum = RunBatch(scoreBuffer, alignmentBuffer);
		scoreBuffer->flush();

		float elapsed = tmr.ET();
		Log.Debug(LOG_CS_DETAILS, "CS Thread %i finished batch (len %i) with %i overflows, length %d (elapsed: %.3fs)", m_TID, m_CurrentBatch.size(), m_Overflows, c_SrchTableBitLen, elapsed);

		NGM.Stats->readsPerSecond = (NGM.Stats->csTime + 1.0f / (elapsed / m_CurrentBatch.size())) / 2.0f;

		NGM.Stats->alignTime = std::max(0.0f, alignmentBuffer->getTime());
		NGM.Stats->scoreTime = std::max(0.0f, scoreBuffer->getTime() - NGM.Stats->alignTime);
		NGM.Stats->csTime = std::max(0.0f, elapsed - NGM.Stats->scoreTime - NGM.Stats->alignTime);

		NGM.Stats->csLength = c_SrchTableBitLen;
		NGM.Stats->csOverflows = m_Overflows;
		NGM.Stats->avgnCRMS = nCRMsSum / m_CurrentBatch.size();

		if (!Config.Exists("search_table_length")) {
			if (m_Overflows <= 5 && !up && c_SrchTableBitLen > 8) {
				SetSearchTableBitLen( c_SrchTableBitLen - 1 );
				Log.Debug(LOG_CS_DETAILS, "Overflow: Switching to %d bits (%d, %d)", c_SrchTableBitLen, c_BitShift, c_SrchTableLen);
			} else if (m_Overflows > m_BatchSize * 0.01f) {
				SetSearchTableBitLen( c_SrchTableBitLen + 1 );
				up = true;
				Log.Debug(LOG_CS_DETAILS, "Overflow: Switching to %d bits (%d, %d)", c_SrchTableBitLen, c_BitShift, c_SrchTableLen);
			}
		}
		m_Overflows = 0;
	}

	delete scoreBuffer;
	scoreBuffer = 0;
	delete alignmentBuffer;
	alignmentBuffer = 0;
	NGM.DeleteAlignment(oclAligner);
	Log.Debug(LOG_CS_DETAILS, "CS Thread %i finished (%i reads processed, %i reads written, %i reads discarded)", m_TID, m_ProcessedReads, m_WrittenReads, m_DiscardedReads);
}
void CS::Init() {
	prefixBasecount = Config.GetInt("kmer", 4, 32);
	prefixBits = prefixBasecount * 2;
	prefixMask = ((ulong) 1 << prefixBits) - 1;

	bool m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);
	if (m_EnableBS) {
		static int const cMutationLocLimit =
		Config.Exists("bs_cutoff") ? Config.GetInt("bs_cutoff") : 6;
		Log.Message("BS mapping enabled. Max. number of A/T per k-mer set to %d", cMutationLocLimit);
	}
}

CS::CS(bool useBuffer) :
		m_CSThreadID((useBuffer) ? (AtomicInc(&s_ThreadCount) - 1) : -1), m_BatchSize(
				cBatchSize / Config.GetInt("qry_avg_len")), m_ProcessedReads(0), m_WrittenReads(
				0), m_DiscardedReads(0), m_EnableBS(false), m_Overflows(0), m_entry(
				0), m_PrefixBaseSkip(0), m_Fallback(false) // cTableLen <= 0 means always use fallback
{

	SetSearchTableBitLen(
	Config.Exists("search_table_length") ?
	Config.GetInt("search_table_length") :
											16);

	m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

	if (m_EnableBS) {
		m_PrefixBaseSkip = (Config.Exists("kmer_skip")) ?
		Config.GetInt("kmer_skip", 0, -1) :
															cPrefixBaseSkip;
	}

	rTable = 0;

	Log.Debug(LOG_CS_DETAILS, "SearchTabLen: %d (%d)", c_SrchTableLen, c_SrchTableBitLen);

#ifdef _DEBUGCS
	ofp2 = fopen("cs-results.txt", "w");
#endif

	currentState = 0;
	tmpSize = 10000;
	tmp = new LocationScore[tmpSize];
}

void CS::AllocRefEntryChain() {
	m_entryCount = m_RefProvider->GetRefEntryChainLength();
	m_entry = new RefEntry[m_entryCount];
}

CS::~CS() {
#ifdef _DEBUGCS
	fclose(ofp2);
#endif
	if (tmp != 0) {
		delete[] tmp;
		tmp = 0;
	}
	if (rTable != 0)
		delete[] rTable;

	if (rList != 0) {
		delete[] rList;
	}

	if (m_entry != 0) {
		delete[] m_entry;
		m_entry = 0;
	}
}

void CS::CheckFallback() {

}

