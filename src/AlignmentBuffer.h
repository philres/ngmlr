#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "GenericReadWriter.h"

#include <map>
#include <list>
#include <algorithm>

#include "FFormatWriter.h"
#include "SAMWriter.h"
#include "BAMWriter.h"

#undef module_name
#define module_name "OUT"

typedef void (*pfDelete)(void*);

char Strand(MappedRead * read);

char const * const cFormatNames[3] = { "Plain Text", "SAM", "BAM" };

class AlignmentBuffer {

private:

	char const * const output_name;
	int const outputformat;
	int const alignmode;
	bool m_EnableBS;
	int const batchSize;
	int const corridor;
	int const refMaxLen;
	MappedRead ** reads;
	int nReads;
	char const * * qryBuffer;
	char const * * refBuffer;
	char * m_DirBuffer;
	int dbLen;
	Align * alignBuffer;
	char * dBuffer;
	char * dummy;

	long tCount;
	long tSum;
	long brokenPairs;
	NGMMutex m_Mutex;

	const char* const m_Filename;
	GenericReadWriter* m_Writer;
	std::map<int, std::list<MappedRead*> > m_EqualScoringBuffer;
	void saveEqualScoring(int id);

	static bool first;

	static bool EqualScoringSortPredicate(MappedRead * lhs, MappedRead * rhs) {
		return (lhs->Strand == rhs->Strand) ? lhs->TLS()->Location.m_Location < rhs->TLS()->Location.m_Location : lhs->Strand < rhs->Strand;
	}

	static bool EqualScoringEquivalentPredicate(MappedRead * lhs, MappedRead * rhs) {
		return (lhs->TLS()->Location.m_RefId == rhs->TLS()->Location.m_RefId)
		&& (lhs->TLS()->Location.m_Location == rhs->TLS()->Location.m_Location);
	}

public:

	static ulong alignmentCount;

	int * RefStartPos;

	IAlignment * aligner;

	AlignmentBuffer(const char* const filename, IAlignment * mAligner) :
			m_Filename(filename),
			batchSize(mAligner->GetAlignBatchSize() / 2),
			output_name(Config.GetString("output")),
			outputformat(NGM.GetOutputFormat()),
			alignmode(Config.GetInt("mode", 0, 1)),
			corridor(Config.GetInt("corridor")),
			refMaxLen(((Config.GetInt("qry_max_len") + corridor) | 1) + 1), aligner(mAligner) {
		RefStartPos = 0;
		tSum = 0;
		tCount = 0;
		brokenPairs = 0;
		m_Writer = 0;
		nReads = 0;


		m_EnableBS = false;
		m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

		int const outputformat = NGM.GetOutputFormat();
//			char const * const output_name = Config.GetString("output");
//		Log.Message("Writing output to %s (Format: %s)", output_name, cFormatNames[outputformat]);
//


						m_Writer =
						(outputformat == 0) ?
						(GenericReadWriter*) new FFormatWriter(
								NGM.getWriter()) :
						(outputformat == 1) ?
						(GenericReadWriter*) new SAMWriter(NGM.getWriter()) :
						(GenericReadWriter*) new BAMWriter(NGM.getWriter());


						///NGMLock(&m_Mutex);
						if(first) {
							Log.Green("Wrinting Prolog!");
							m_Writer->WriteProlog();
							first = false;
						}
						///NGMUnlock(&m_Mutex);

//		int const batchSize = NGM.Aligner()->GetAlignBatchSize();
		//m_Writer = NGM.getWriter(output_name);


//		NGMInitMutex(&m_Mutex);

		Log.Verbose("Alignment batchsize = %i", batchSize);

		reads = new MappedRead*[batchSize];

		qryBuffer = new char const *[batchSize];
		refBuffer = new char const *[batchSize];

		for (int i = 0; i < batchSize; ++i) {
			refBuffer[i] = new char[refMaxLen];
		}

		m_DirBuffer = new char[batchSize];

		alignBuffer = new Align[batchSize];
		dbLen = std::max(1, Config.GetInt("qry_max_len")) * 8;
		dBuffer = new char[dbLen];

		dummy = new char[refMaxLen];
		memset(dummy, '\0', refMaxLen);
		//dummy[Config.GetInt("qry_max_len") - 1] = '\0';

	}

	virtual ~AlignmentBuffer() {
		delete m_Writer;
		Log.Message("Freeing resources...");
		delete[] m_DirBuffer;
		m_DirBuffer = 0;

		delete[] dummy;
		dummy = 0;

		for (int i = 0; i < batchSize; ++i) {
			delete[] refBuffer[i];
			refBuffer[i] = 0;
		}
		delete[] qryBuffer;
		delete[] refBuffer;
		delete[] alignBuffer;

		delete[] reads;
		delete[] dBuffer;

		Log.Message("Finalizing output...");
		//m_Writer->WriteEpilog();

		//delete m_Writer;

		Log.Message("Output finished");
	}

	void DoRun();


	int GetStage() const {
		return 4;
	}

	inline const char* GetName() const {
		return "Output";
	}

	void addRead(MappedRead * read);
	void flush();


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

		//NGMLock(&m_Mutex);
//		NGM.AquireOutputLock();
//#pragma omp critical
	//	{
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
//		NGM.ReleaseOutputLock();
		//NGMUnlock(&m_Mutex);
	//}

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
};

#endif
