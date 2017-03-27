/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "CS.h"

#include <memory.h>
#include <stdlib.h>
#include <cmath>

#include "NGM.h"
#include "Timing.h"
#include "AlignmentBuffer.h"

#undef module_name
#define module_name "CS"

bool stdoutPrintRunningTime = false;

//static const int cInvalidLocation = -99999999;
static const SequenceLocation sInvalidLocation(9223372036854775808u, 0, false);

volatile int CS::s_ThreadCount = 0;

bool up = false;
int kCount = 0;

//Default values
int x_SrchTableBitLen = 24;
//Number of bases per batch
//Default batch size
//int const cBatchSize = 1800000;
// Reduced batch size (for low read number PacBio samples)
int const cBatchSize = 10;
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

	prefixBasecount = Config.getKmerLength();
	prefixBits = prefixBasecount * 2;
	prefixMask = ((ulong) 1 << prefixBits) - 1;
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

		for (int i = 0; i < n; ++i) {
			uloc loc = cur->getRealLocation(cur->ref[i]);
			cs->AddLocationStd(GetBin(loc - correction), cur->reverse, weight);
			//cs->AddLocationStd(GetBin(loc - (int) (correction * 0.95f)),
			//		cur->reverse, weight);
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
	ss << std::string(Config.getOutputFile()) << "_kmer-profile/" << read->ReadId
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

	static float const mink = Config.getMinKmerHits();

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
			toInsert->Location.m_Location = ResolveBin(temp.m_Location);
			toInsert->Location.setReverse(false);
		}
		if (temp.rScore >= mi_Threshhold) {
			LocationScore * toInsert = &tmp[index++];
			toInsert->Score.f = temp.rScore;
			toInsert->Location.m_Location = ResolveBin(temp.m_Location);
			toInsert->Location.setReverse(true);
		}
	}
	static int const maxScores = Config.getMaxCMRs();
	if (index < maxScores)
		read->AllocScores(tmp, index);

	return index;
}

int CS::CollectResultsFallback(MappedRead * read) {
	throw "Not implemented (CollectResultsFallback)";
	return 0;
}

void CS::SendToBuffer(MappedRead * read, ScoreBuffer * sw,
		AlignmentBuffer * out) {

	if (read == 0)
		return;

	int count = read->numScores();
	if (read->group != 0) {
		/*
		 * Read is a long read (long enough for splitting)
		 */
		if (count == 0) {
			/*
			 * Current sub-read does not have mapping locations
			 */
			read->Calculated = 0;
			read->mappingQlty = 0;
			read->group->readsFinished += 1;

			if (read->group->readsFinished == read->group->readNumber) {
				out->processLongReadLIS(read->group);
			}
		} else {
			/*
			 * Current sub-read has mapping locations
			 *  -> compute scores
			 */
			read->Calculated = 0;
			sw->addRead(read, count);
		}

	} else {
		/*
		 * Read length is too short for splitting
		 */
		if (count == 0) {
			out->WriteRead(read, false);
		} else {
			read->Calculated = 0;
			sw->scoreShortRead(read);
		}
	}
}

void CS::Cleanup() {
	NGM.ReleaseWriter();
}

int CS::RunRead(MappedRead * currentRead, PrefixIterationFn pFunc,
		ScoreBuffer * sw, AlignmentBuffer * out) {

	int nScoresSum = 0;

	m_CurrentSeq = currentRead->ReadId;
	//m_CurrentSeq = i;

	//currentState = 2 * i;
	currentState++;
	rListLength = 0;
	maxHitNumber = 0.0f;
	currentThresh = 0.0f;
	//weightSum = 0.0f;
	kCount = 0;

	m_CurrentReadLength = currentRead->length;

	++m_ProcessedReads;
	bool fallback = m_Fallback;

	char const * const qrySeq = currentRead->Seq;
	uloc qryLen = currentRead->length;

	//		SetSearchTableBitLen(24);

	if (!fallback) {
		try {
			hpoc = c_SrchTableLen * 0.333f;
			PrefixIteration(qrySeq, qryLen, pFunc, 0, 0, this,
					m_PrefixBaseSkip);
			int nScores = CollectResultsStd(currentRead);
			nScoresSum += nScores;
		} catch (int overflow) {
			++m_Overflows;

			// Initiate fallback
			fallback = true;
		}
	}
	if (fallback) {
		int c_SrchTableBitLenBackup = c_SrchTableBitLen;
		int x = 2;
		while (fallback && (c_SrchTableBitLenBackup + x) <= 20) {
			fallback = false;
			try {

				SetSearchTableBitLen(c_SrchTableBitLenBackup + x);

				rListLength = 0;
				//currentState = 2 * i + 1;
				currentState++;
				maxHitNumber = 0.0f;
				currentThresh = 0.0f;

				hpoc = c_SrchTableLen * 0.777f;
				PrefixIteration(qrySeq, qryLen, pFunc, 0, 0,
						this, m_PrefixBaseSkip);
				nScoresSum += CollectResultsStd(currentRead);

			} catch (int overflow) {
				//					Log.Message("Overflow in fallback for read %s!", currentRead->name);
				fallback = true;
				x += 1;
			}
		}
		if (fallback) {
			Log.Message("Couldn't find candidate for read %s (too many candidates)", currentRead->name);
		}
		SetSearchTableBitLen(c_SrchTableBitLenBackup);
	}

	SendToBuffer(currentRead, sw, out);
	return nScoresSum;
}

int CS::RunBatch(ScoreBuffer * sw, AlignmentBuffer * out) {
	PrefixIterationFn pFunc = &CS::PrefixSearch;

	int nScoresSum = 0;
	for (size_t i = 0; i < m_CurrentBatch.size(); ++i) {

		MappedRead * currentRead = m_CurrentBatch[i];
		nScoresSum += RunRead(currentRead, pFunc, sw, out);
	}
	return nScoresSum;
}

void CS::DoRun() {


	NGM.AquireOutputLock();
	oclAligner = NGM.CreateAlignment(0);

	alignmentBuffer = new AlignmentBuffer(Config.getOutputFile());
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
	m_CsSensitivity = Config.getSensitivity();

	m_RefProvider = 0;

	int batchId = 0;

	Timer tmr;
	tmr.ST();
	m_RefProvider = NGM.GetRefProvider(m_TID);
	AllocRefEntryChain();
	while (NGM.ThreadActive(m_TID, GetStage()) && ((m_CurrentBatch = NGM.GetNextReadBatch(m_BatchSize)), (m_CurrentBatch.size() > 0))) {
		Log.Debug(LOG_CS_DETAILS, "CS Thread %i got batch (len %i)", m_TID, m_CurrentBatch.size());
//		Log.Message("CS Thread %i got batch %d (len %i)", m_TID, ++batchId, m_CurrentBatch.size());

		int nCRMsSum = RunBatch(scoreBuffer, alignmentBuffer);
		Log.Verbose("CS Thread %i flushing scoreBuffer for batch %d", m_TID, batchId);
		scoreBuffer->flush();

		float elapsed = tmr.ET();
		Log.Debug(LOG_CS_DETAILS, "CS Thread %i finished batch (len %i) with %i overflows, length %d (elapsed: %.3fs)", m_TID, m_CurrentBatch.size(), m_Overflows, c_SrchTableBitLen, elapsed);
//		Log.Message("CS Thread %i finished batch %d (len %i) with %i overflows, length %d (elapsed: %.3fs)", m_TID, batchId, m_CurrentBatch.size(), m_Overflows, c_SrchTableBitLen, elapsed);

		float alignmentTime = alignmentBuffer->getAlignTime();
		float longReadProcessingTime = alignmentBuffer->getProcessTime() - alignmentTime;
		float shortReadScoringTime = scoreBuffer->getTime() - alignmentBuffer->getProcessTime();
		float kmerSearchTime = elapsed - scoreBuffer->getTime();

//		if(stdoutPrintRunningTime) {
//			fprintf(stdout, "%f\t%f\t%f\t%f\t%f\n", elapsed,
//					scoreBuffer->getTime(),
//					alignmentBuffer->getTime(),
//					alignmentBuffer->getProcessTime(),
//					alignmentBuffer->getAlignTime());
//			fprintf(stdout, "%f (%d)\t%f(%d)\t%f(%d)\t%f(%d)\n",
//					kmerSearchTime, (int)(kmerSearchTime / elapsed * 100.0f),
//					shortReadScoringTime, (int)(shortReadScoringTime / elapsed * 100.0f),
//					longReadProcessingTime, (int)(longReadProcessingTime / elapsed * 100.0f),
//					alignmentTime, (int)(alignmentTime / elapsed * 100.0f));
//			fflush(stdout);
//		}

		NGM.Stats->alignTime = (NGM.Stats->alignTime + (int)(alignmentTime / elapsed * 100.0f)) / 2;
		NGM.Stats->scoreTime = (NGM.Stats->scoreTime + (int)(shortReadScoringTime / elapsed * 100.0f)) / 2;
		NGM.Stats->csTime = (NGM.Stats->csTime + (int)(kmerSearchTime / elapsed * 100.0f)) / 2;

		NGM.Stats->csLength = c_SrchTableBitLen;
		NGM.Stats->csOverflows = m_Overflows;
		NGM.Stats->avgnCRMS = nCRMsSum / m_CurrentBatch.size();

		if (Config.getCsSearchTableLength() == 0) {
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
//	Log.Message("CS Thread %i finished (%i reads processed, %i reads written, %i reads discarded)", m_TID, m_ProcessedReads, m_WrittenReads, m_DiscardedReads);
}
void CS::Init() {
	prefixBasecount = Config.getKmerLength();
	prefixBits = prefixBasecount * 2;
	prefixMask = ((ulong) 1 << prefixBits) - 1;
}

CS::CS(bool useBuffer) :
		m_CSThreadID((useBuffer) ? (AtomicInc(&s_ThreadCount) - 1) : -1), m_BatchSize(
				cBatchSize), m_ProcessedReads(0), m_DiscardedReads(
				0), m_Overflows(0), m_entry(0), m_PrefixBaseSkip(
				0), m_Fallback(false), alignmentBuffer(0), c_SrchTableLen(0), rListLength(0), m_CurrentSeq(0), c_SrchTableBitLen(0), oclAligner(0),
				m_CurrentReadLength(0), hpoc(0), m_CsSensitivity(0.0f), maxHitNumber(0.0f), rList(0), m_entryCount(0), c_BitShift(0), m_RefProvider(0), currentThresh(0.0f)
{

	SetSearchTableBitLen(Config.getCsSearchTableLength() != 0 ? Config.getCsSearchTableLength() : 16);

	rTable = 0;

	Log.Debug(LOG_CS_DETAILS, "SearchTabLen: %d (%d)", c_SrchTableLen, c_SrchTableBitLen);

	currentState = 0;
	tmpSize = 10000;
	tmp = new LocationScore[tmpSize];
}

void CS::AllocRefEntryChain() {
	m_entryCount = m_RefProvider->GetRefEntryChainLength();
	m_entry = new RefEntry[m_entryCount];
}

CS::~CS() {
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

