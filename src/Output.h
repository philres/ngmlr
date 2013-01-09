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

	Output(const char* const filename) :
			m_Filename(filename) {
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
			//Correct position
			int refCount = SequenceProvider.GetRefCount();
			SequenceLocation loc = read->TLS()->Location;
			int i = 0;
			while (i < refCount && loc.m_Location >= SequenceProvider.GetRefStart(i)) {
				i += (NGM.DualStrand()) ? 2 : 1;
			}
			if (i == 0) {
				Log.Message("%s (%d) - %s: %s, %s", read->name, read->ReadId, read->Seq, read->Buffer1, read->Buffer2);
				Log.Error("Couldn't resolve mapping position: %u!", loc.m_Location);
				Fatal();
			}
			i -= (NGM.DualStrand()) ? 2 : 1;

			loc.m_Location -= SequenceProvider.GetRefStart(i);
			loc.m_RefId = i;
			read->TLS()->Location = loc;

			if (read->Strand == '-') {
				if (read->qlty != 0)
					std::reverse(read->qlty, read->qlty + strlen(read->qlty));
			}
		}
		NGM.AquireOutputLock();
		if (NGM.Paired() && read->Paired != 0) {
			if (read->Paired->HasFlag(NGMNames::DeletionPending)) {
				m_Writer->WritePair(read, read->Paired);
				//NGM.AddWrittenRead(read->ReadId);
				//NGM.AddWrittenRead(read->Paired->ReadId);
			}
		} else {
			m_Writer->WriteRead(read, mapped);
			//NGM.AddWrittenRead(read->ReadId);
		}
		NGM.ReleaseOutputLock();
	}

private:
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
