#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "NGMTask.h"
#include "SW.h"
#include "GenericReadWriter.h"

#include <map>
#include <list>
#include <algorithm>

#undef module_name
#define module_name "OUT"

typedef void (*pfDelete)(void*);

char Strand(MappedRead * read);

class Output: public NGMTask {
public:

	static ulong alignmentCount;

	int * RefStartPos;

	Output(const char* const filename) :
			m_Filename(filename) {
		RefStartPos = 0;
		tSum = 0;
		tCount = 0;
		brokenPairs = 0;
	}

	void DoRun();

	int GetStage() const {
		return 4;
	}

	inline const char* GetName() const {
		return "Output";
	}

	void SaveRead(MappedRead* read, bool mapped = true) {
		if (mapped) {

			static int refCount = SequenceProvider.GetRefCount();
			if (RefStartPos == 0) {
				RefStartPos = new int[refCount / ((NGM.DualStrand()) ? 2 : 1)];
				int i = 0;
				int j = 0;
				while (i < refCount/* && loc.m_Location >= SequenceProvider.GetRefStart(i)*/) {
					RefStartPos[j++] = SequenceProvider.GetRefStart(i);
					i += (NGM.DualStrand()) ? 2 : 1;
				}
			}

			//Correct position
			SequenceLocation loc = read->TLS()->Location;
			int * upper = std::upper_bound(RefStartPos, RefStartPos + (refCount / ((NGM.DualStrand()) ? 2 : 1)), loc.m_Location);
			std::ptrdiff_t refId = ((upper - 1) - RefStartPos) * ((NGM.DualStrand()) ? 2 : 1);
			loc.m_Location -= *(upper - 1);
			loc.m_RefId = refId;
			read->TLS()->Location = loc;

			if (read->Strand == '-') {
				if (read->qlty != 0)
				std::reverse(read->qlty, read->qlty + strlen(read->qlty));
			}
		}
		NGM.AquireOutputLock();
		if (NGM.Paired() && read->Paired != 0) {
			if (read->Paired->HasFlag(NGMNames::DeletionPending)) {

				if(read->hasCandidates() && read->Paired->hasCandidates()) {
					LocationScore * ls1 = read->TLS();
					LocationScore * ls2 = read->Paired->TLS();
					int distance =
					(ls2->Location.m_Location > ls1->Location.m_Location) ?
					ls2->Location.m_Location - ls1->Location.m_Location + ls2->Read->length :
					ls1->Location.m_Location - ls2->Location.m_Location + ls1->Read->length;

					//int distance = abs(read->TLS()->Location.m_Location - read->Paired->TLS()->Location.m_Location);

					tCount += 1;
					if(ls1->Location.m_RefId != ls2->Location.m_RefId || distance < _NGM::sPairMinDistance || distance > _NGM::sPairMaxDistance || read->Strand == read->Paired->Strand) {
						read->SetFlag(NGMNames::PairedFail);
						read->Paired->SetFlag(NGMNames::PairedFail);
						brokenPairs += 1;
					} else {
						//Log.Message("%d", distance);
						tSum += distance;
					}
				}
				m_Writer->WritePair(read, read->Paired);
			}
		} else {
			m_Writer->WriteRead(read, mapped);
		}
		NGM.ReleaseOutputLock();

		if(tCount % 1000 == 0) {
//			Log.Message("%d %d %d %f", tSum, tCount, brokenPairs, tSum * 1.0f / (tCount));
			NGM.Stats->validPairs = (tCount - brokenPairs) * 100.0f / tCount;
			NGM.Stats->insertSize = tSum * 1.0f / (tCount - brokenPairs);
//			Log.Message("%f %f", NGM.Stats->validPairs, NGM.Stats->insertSize);
		} else {
			//tSum = 0;
			//tCount = 1;
			//brokenPairs = 0;
		}
	}

private:

	long tCount;
	long tSum;
	long brokenPairs;

	const char* const m_Filename;
	GenericReadWriter* m_Writer;
	std::map<int, std::list<MappedRead*> > m_EqualScoringBuffer;
	void saveEqualScoring(int id);

	static bool EqualScoringSortPredicate(MappedRead * lhs, MappedRead * rhs) {
		return (lhs->Strand == rhs->Strand) ? lhs->TLS()->Location.m_Location < rhs->TLS()->Location.m_Location : lhs->Strand < rhs->Strand;
	}

	static bool EqualScoringEquivalentPredicate(MappedRead * lhs, MappedRead * rhs) {
		return (lhs->TLS()->Location.m_RefId == rhs->TLS()->Location.m_RefId)
		&& (lhs->TLS()->Location.m_Location == rhs->TLS()->Location.m_Location);
	}
}
;

#endif
