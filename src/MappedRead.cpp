#include "MappedRead.h"
#include "Log.h"
#include "Config.h"

#undef module_name
#define module_name "READ"

#ifdef INSTANCE_COUNTING
volatile int MappedRead::sInstanceCount = 0;
volatile int LocationScore::sInstanceCount = 0;
#endif

volatile int MappedRead::maxSeqCount = 0;

MappedRead::MappedRead(int const readid) :
		ReadId(readid), Calculated(-1), TopScore(-1), Lock(0), EqualScoringID(0), EqualScoringCount(1), Scores(), nScores(0), iScores(0), Paired(
				0), Status(0), Identity(0), NM(-1), Strand('?'), QStart(0), QEnd(0), mappingQlty(255), RevSeq(0), Seq(0), qlty(0), name(0), Buffer1(
				0), Buffer2(0) {
#ifdef INSTANCE_COUNTING
	AtomicInc(&sInstanceCount);
	maxSeqCount = std::max(sInstanceCount, maxSeqCount);
#endif
}

void MappedRead::AllocBuffers() {
	static int const qryMaxLen = Config.GetInt("qry_max_len");
	Buffer1 = new char[std::max(1, qryMaxLen) * 4];
	Buffer2 = new char[std::max(1, qryMaxLen) * 4];
	*(int*) Buffer1 = 0x212121;
	*(int*) Buffer2 = 0x212121;
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

char const * MappedRead::computeReverseSeq() {
	if (RevSeq == 0) {
		static int const qryMaxLen = Config.GetInt("qry_max_len");
		RevSeq = new char[qryMaxLen];
		memset(RevSeq, 0, qryMaxLen);

		char * fwd = Seq;
		char * rev = RevSeq + length - 1;

		for (int i = 0; i < length; ++i) {
			*rev-- = cpl(*fwd++);
			//Seq[i] = cpl(RevSeq[length - i - 1]);
		}
	}
	return RevSeq;
}

void MappedRead::DeleteReadSeq() {
	if (Seq != 0)
		delete[] Seq;
	Seq = 0;
	if (RevSeq != 0)
		delete[] RevSeq;
	RevSeq = 0;
}

MappedRead::~MappedRead() {
	//for (size_t i = 0; i < Scores.size(); ++i)
	//	delete Scores[i];
	if (Scores != 0) {
		delete[] Scores;
		Scores = 0;
		iScores = 0;
		nScores = 0;
	}
	if (Buffer1 != 0)
		delete[] Buffer1;
	if (Buffer2 != 0)
		delete[] Buffer2;
	DeleteReadSeq();
	if (qlty != 0) {
		delete[] qlty;
		qlty = 0;
	}
	if (name != 0)
		delete[] name;
	name = 0;

#ifdef INSTANCE_COUNTING
	AtomicDec(&sInstanceCount);
#endif
}

void MappedRead::AllocScores(LocationScore * tmp, int const n) {
	nScores = n;
	iScores = n;
//	Log.Message("Allocating %d scores.", nScores);
	if (nScores > 0) {
		Scores = new LocationScore[nScores];
		memcpy(Scores, tmp, n * sizeof(LocationScore));
	}
}

//LocationScore * MappedRead::AddScore(LocationScore const & score) {
//	LocationScore * toInsert = new LocationScore(score);
//	toInsert->Read = this;
//	Scores.push_back(toInsert);
//	return toInsert;
//}

LocationScore * MappedRead::AddScore(float const score, uint const loc, bool const reverse) {
	Log.Error("AddScore");
	throw "";
	iScores += 1;
	LocationScore * toInsert = &Scores[iScores - 1];
	toInsert->Score.f = score;
	toInsert->Location.m_Location = loc;
	toInsert->Location.m_RefId = reverse;
	toInsert->Read = this;
	return toInsert;
	//LocationScore * toInsert = new LocationScore(loc, reverse, score, this);
	//toInsert->Read = this;
	//Scores.push_back(toInsert);
	//return toInsert;
}

void MappedRead::clearScores(bool const keepTop) {
	if (TopScore != -1 && Scores != 0) {
		LocationScore tmp = *TLS();

		delete[] Scores;

		if (keepTop) {
			Scores = new LocationScore[1];
			Scores[0] = tmp;
			TopScore = 0;
			iScores = 1;
			nScores = 1;
		} else {
			TopScore = -1;
			iScores = 0;
			nScores = 0;
			Scores = 0;
		}
	}

}

LocationScore * MappedRead::TLS() const {
	if (TopScore == -1) {
		Log.Error("Top score accessed without calculation (Read %i, ESID %i, ESC %i)", ReadId, EqualScoringID, EqualScoringCount);
		throw "";
		return 0;
	}
	return &Scores[TopScore];
}

int MappedRead::numScores() const {
	return iScores;
}

bool MappedRead::hasCandidates() const {
	return iScores > 0;
}

void PrintScores(MappedRead const * read) {
//	Log.Green("Scores for read %i", read->ReadId);
//	for (size_t i = 0; i < read->Scores.size(); ++i) {
//		if (read->Scores[i] == 0)
//			Log.Message("<0>");
//		else
//			Log.Message("Score %.2f for <%i, %i>", read->Scores[i]->Score.f, read->Scores[i]->Location.m_Location, read->Scores[i]->Location.m_RefId);
//	}
}

//void MappedRead::TopN(uint n) {
//	if (Scores.size() == 0)
//		return;
//
//	if (Calculated > 0 && Calculated != (int) Scores.size())
//		Log.Error("TopN called for inconsistent scores may lead to unexpected results.");
//
//	std::sort(Scores.begin(), Scores.end(), SortPred); // Sort scores from best to worst
//	Unique(); // Remove double scores
//	Scores.erase(std::remove_if(Scores.begin(), Scores.end(), IsZero), Scores.end()); // remove 0-pointer from last step
//
//	// remove from back all but top n
//	while (Scores.size() > n) {
//		delete Scores[Scores.size() - 1];
//		Scores.pop_back();
//	}
//
//	// realloc vector
//	std::vector<LocationScore*>(Scores).swap(Scores);
//
//	TopScore = 0;
//	Calculated = n;
//}

//void MappedRead::Unique() {
//	uint ci = 0;
//	for (uint i = 1; i < Scores.size(); ++i) {
//		if (UniquePred(Scores[ci], Scores[i])) {
//			delete Scores[i];
//			Scores[i] = 0;
//		} else
//			ci = i;
//	}
//}
