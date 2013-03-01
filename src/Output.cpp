#include <stdio.h>

#include "Timing.h"
#include "Output.h"
#include "Debug.h"

#include <string.h>

ulong Output::alignmentCount = 0;

//static inline void rev(char * s);

void Output::saveEqualScoring(int id) {
	typedef std::list<MappedRead*> TReadList;
	TReadList reads = m_EqualScoringBuffer[id];
	//							if (Config.GetInt("bowtie_mode") == 0) {
	reads.sort(EqualScoringSortPredicate);
	//reads.unique(EqualScoringEquivalentPredicate);
	int esc = reads.size();
	TReadList::iterator itr = reads.begin();
	TReadList::iterator citr = itr;
	++itr;
	while (itr != reads.end()) {
		if (EqualScoringEquivalentPredicate(*citr, *itr)) {
//									if (esc == 1) {
//										(**itr).DeleteReadSeq();
//									}
			//delete *itr;
			NGM.GetReadProvider()->DisposeRead(*itr);
			*itr = 0;
			--esc;
		} else
		citr = itr;
		++itr;
	}
	int esid = 0;
	for (TReadList::iterator itr = reads.begin(); itr != reads.end(); ++itr) {
		if (*itr != 0) {
			(**itr).EqualScoringID = esid++;
			(**itr).EqualScoringCount = esc;
			SaveRead(*itr);
//									if (esid == cur_read->EqualScoringCount) {
//										(**itr).DeleteReadSeq();
//									}
			NGM.GetReadProvider()->DisposeRead(*itr);
			//delete *itr;

			*itr = 0;
		}
	}
	//std::for_each(reads.begin(), reads.end(), std::bind1st( std::mem_fun(&Output::SaveRead), this ) );
	//							} else {
	//								SaveRead(cur_read);
	//								for (TReadList::iterator itr = reads.begin(); itr != reads.end(); ++itr) {
	//									NGM.GetReadProvider()->DisposeRead(*itr);
	//								}
	//							}
	reads.clear();
	m_EqualScoringBuffer[id].clear();
	std::map<int, std::list<MappedRead*> >::iterator it;
	it = m_EqualScoringBuffer.find(id);
	m_EqualScoringBuffer.erase(it);
}

void Output::flush() {
	DoRun();
	nReads = 0;
}

void Output::addRead(MappedRead * read) {
	if (!read->hasCandidates()) {
		SaveRead(read, false);
	} else {
		reads[nReads++] = read;
		if (nReads == 4096) {
			DoRun();
			nReads = 0;
		}
	}
}

void Output::DoRun() {

	int count = nReads;
	Timer tmr;
	tmr.ST();

	if (count > 0) {
		alignmentCount += count;
		for (int i = 0; i < count; ++i) {
			MappedRead * cur_read = reads[i];
			if (cur_read->hasCandidates()) {
				cur_read->Strand = Strand(cur_read);

				// Bei Mapping am Minus-Strang, Position bez�glich +-Strang reporten
				if (cur_read->Strand == '-') {
					qryBuffer[i] = cur_read->RevSeq;
					//Log.Message("Rev Seq: %s", qryBuffer[i]);

					// RefId auf +-Strang setzen
					--cur_read->TLS()->Location.m_RefId;
					if (NGM.Paired()) {
						m_DirBuffer[i] = !(cur_read->ReadId & 1);
					} else {
						m_DirBuffer[i] = 1;
					}

				} else {
					qryBuffer[i] = cur_read->Seq;
					if (NGM.Paired()) {
						m_DirBuffer[i] = cur_read->ReadId & 1; //0 if first pair
					} else {
						m_DirBuffer[i] = 0;
					}
				}

				SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), cur_read->TLS()->Location.m_RefId,
						cur_read->TLS()->Location.m_Location - (corridor >> 1), refMaxLen);

				cur_read->AllocBuffers();
				alignBuffer[i].pBuffer1 = cur_read->Buffer1;
				alignBuffer[i].pBuffer2 = cur_read->Buffer2;

			} else {
				Log.Warning("Unmapped read submitted to alignment computation!");

				Fatal();
				char * refDummy = const_cast<char *>(refBuffer[i]);
				memset(refDummy, '\0', refMaxLen);
				qryBuffer[i] = dummy;//SequenceProvider.GetQrySequence(cur_read->ReadId);

				alignBuffer[i].pBuffer1 = dBuffer;
				alignBuffer[i].pBuffer2 = dBuffer + dbLen / 2;

			}
			//Log.Message("%s: %s", cur_read->name, refBuffer[i]);
			//Log.Message("%s: %s", cur_read->name, qryBuffer[i]);
		}

		Log.Verbose("Thread %i invoking alignment (count = %i)", 0, count);
		Timer x;
		x.ST();
		int aligned =
				aligner->BatchAlign(alignmode | (std::max(outputformat, 1) << 8), count, refBuffer, qryBuffer, alignBuffer,
						(m_EnableBS) ? m_DirBuffer : 0);
		//int aligned = count;

		if (aligned == count) {
			Log.Verbose("Output Thread %i finished batch (Size = %i, Elapsed: %.2fs)", 0, count, x.ET());
//				Log.Message("Align done");
		} else {
			Log.Error("Error aligning outputs (%i of %i aligned)", aligned, count);
		}

		for (int i = 0; i < aligned; ++i) {
			MappedRead * cur_read = reads[i];

			Log.Verbose("Process aligned read %i,%i (%s)", i, cur_read->ReadId, cur_read->name);
			int id = cur_read->ReadId;

			if (cur_read->hasCandidates()) {
				cur_read->TLS()->Location.m_Location +=
						alignBuffer[i].PositionOffset - (corridor >> 1);
				//cur_read->TLS()->Score.f = alignBuffer[i].Score; // TODO: Align liefert keine Scores
				cur_read->Identity = alignBuffer[i].Identity;
				cur_read->NM = alignBuffer[i].NM;
				cur_read->QStart = alignBuffer[i].QStart;
				cur_read->QEnd = alignBuffer[i].QEnd;

//					Log.Message("%s: %d %d %f -> %s, %s", cur_read->name, cur_read->QStart, cur_read->QEnd, cur_read->Identity, cur_read->Buffer1, cur_read->Buffer2);

#ifdef _DEBUGOUT
					Log.Message("Read:   %s", cur_read->name);
					Log.Message("Score: %f", cur_read->TLS()->Score.f);
					Log.Message("Seq:    %s", cur_read->Seq);
					Log.Message("CIGAR:  %s", cur_read->Buffer1);
					Log.Message("MD:     %s", cur_read->Buffer2);
#endif

				//TODO: fix equal scoring
				if (false && cur_read->EqualScoringCount > 1) {
					Log.Error("Should not be here! Equal scoring not supported at the moment.");
					//Reads maps to several locations with an equal score
					m_EqualScoringBuffer[id].push_back(cur_read);
					if (m_EqualScoringBuffer[id].size() == cur_read->EqualScoringCount) {
						//Alignments for all equal scoring positions of this reads are computed
						//Warning: TAKE CARE OF UNMAPPABLE READ (INVALID CIGAR STRING0
						saveEqualScoring(id);
					}
				} else {
					SaveRead(cur_read, alignBuffer[i].Score != -1.0f);
					NGM.GetReadProvider()->DisposeRead(cur_read);
				}
			} else {
				Log.Warning("Unmapped read detected during alignment computation!");
				//throw "geh scheiszen";
				//SaveRead(cur_read, false);
				//NGM.GetReadProvider()->DisposeRead(cur_read);
				//Fatal();
			}
		} // for
		Log.Verbose("Output Thread %i finished batch in %.2fs", 0, tmr.ET());
	} // if count > 0
	else {
		Log.Message("Nothing to do...waiting");
	}
}

//static inline char cpl(char c)
//{
//	if (c == 'A')
//		return 'T';
//	else if (c == 'T')
//		return 'A';
//	else if (c == 'C')
//		return 'G';
//	else if (c == 'G')
//		return 'C';
//	else
//		return c;
//}
//
//// swaps two bases and complements them
//static inline void rc(char & c1, char & c2)
//{
//	char x = c1;
//	c1 = cpl(c2);
//	c2 = cpl(x);
//}
//
//// Reverse-complements a 0-terminated string in place
//static inline void rev(char * s)
//{
//	char * end = s + strlen(s) - 1;
//
//	while (s < end)
//		rc(*s++, *end--);
//
//	if (s == end)
//		*s = cpl(*s);
//}

char Strand(MappedRead * read) {
	static bool sDualStrand = NGM.DualStrand();
	return (sDualStrand && (read->TLS()->Location.m_RefId & 1)) ? '-' : '+';
}
