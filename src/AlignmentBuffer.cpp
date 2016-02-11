#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>
#include <cmath>

#include "seqan/EndToEndAffine.h"
#include "OutputReadBuffer.h"
#include "Timing.h"

ulong AlignmentBuffer::alignmentCount = 0;
bool AlignmentBuffer::first = true;

void AlignmentBuffer::flush() {
	DoRun();
	nReads = 0;
}

void AlignmentBuffer::debugAlgnFinished(MappedRead * read) {
	Log.Debug(32, "READ_%d\tALGN\tAll alignments computed (%d)", read->ReadId, read->numScores());

	if(read->numScores() > 0) {
		for(int i = 0; i < read->numScores(); ++i) {

			LocationScore score = read->Scores[i];
			Align align = read->Alignments[i];

			SequenceLocation loc = score.Location;
			SequenceProvider.convert(loc);

			int refNameLength = 0;
			//TODO_GENOMESIZE: Re-enable me
			//Log.Debug(128, "READ_%d\tALGN_RESULTS\tCMR_%d\t%f\t%f\t%d\t%s\t%s\t%d\t%s", read->ReadId, i, score.Score.f, align.Identity, align.NM, align.pBuffer1, align.pBuffer2, loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), refNameLength));
		}

	}

#ifdef _DEBUGCMRS
	SequenceLocation rloc = SequenceProvider.convert(cur_read, cur_read->Scores[scoreId].Location.m_Location);
	int refNameLength = 0;
	fprintf(cmrBed, "%s\t%d\t%d\t%s_%d\t%f\t%c\n", SequenceProvider.GetRefName(rloc.getrefId(), refNameLength), rloc.m_Location - (corridor >> 1), rloc.m_Location - (corridor >> 1) + refMaxLen, cur_read->name, scoreId, cur_read->Scores[scoreId].Score.f, (rloc.isReverse()) ? '-' : '+');
#endif
}

void AlignmentBuffer::addRead(MappedRead * read, int scoreID) {
	if (argos) {
		SaveRead(read, read->hasCandidates());
	} else {
		if (!read->hasCandidates() || read->mappingQlty < min_mq) {
			//If read has no CMRs or mapping quality is lower than min mapping quality, output unmapped read
			//read->clearScores(-1);
			SaveRead(read, false);
		} else {
			Log.Debug(512, "READ_%d\tALGN_BUFFER\tCMR_%d %f (location %llu) added to alignment buffer at position %d", read->ReadId, scoreID, read->Scores[scoreID].Score.f, read->Scores[scoreID].Location.m_Location, nReads);
			//add alignment computations to buffer. if buffer is full, submit to CPU/GPU
			reads[nReads].scoreId = scoreID;
			reads[nReads++].read = read;
			if (nReads == batchSize) {
				DoRun();
				nReads = 0;
			}
		}
	}
}

void AlignmentBuffer::DoRun() {

	int count = nReads;

	if (count > 0) {
		Log.Debug(32, "INFO\tALGN\tSubmitting %d alignment computations.", count);
		Timer tmr;
		tmr.ST();
		alignmentCount += count;
		for (int i = 0; i < count; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;

			assert(cur_read->hasCandidates());

			//Initialize
			if (cur_read->Scores[scoreID].Location.isReverse()) {
				qryBuffer[i] = cur_read->RevSeq;

				if (cur_read->Paired != 0) {
					m_DirBuffer[i] = !(cur_read->ReadId & 1);
				} else {
					m_DirBuffer[i] = 1;
				}

			} else {
				qryBuffer[i] = cur_read->Seq;
				if (cur_read->Paired != 0) {
					m_DirBuffer[i] = cur_read->ReadId & 1; //0 if first pair
				} else {
					m_DirBuffer[i] = 0;
				}
			}

			//decode reference sequence
			if (!SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
							cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen)) {
//							cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor)) {
				Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
				//Log.Warning("Read sequence: %s", cur_read->Seq);
				memset(const_cast<char *>(refBuffer[i]), 'N', refMaxLen);
			}
//			//decode reference sequence
//			SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
//					cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen);

			//initialize arrays for CIGAR and MD string
			Log.Error("Should not be here. AlignmentBuffer DoRun");
			Fatal();
			static int const qryMaxLen = Config.GetInt("qry_max_len");
			alignBuffer[i].pBuffer1 = new char[std::max(1, qryMaxLen) * 4];
			alignBuffer[i].pBuffer2 = new char[std::max(1, qryMaxLen) * 4];
			*(int*) alignBuffer[i].pBuffer1 = 0x212121;
			*(int*) alignBuffer[i].pBuffer2 = 0x212121;

			//Log.Message("Ref:  %s\nRead: %s", refBuffer[i], qryBuffer[i]);

		}

		//start alignment
		int aligned = aligner->BatchAlign(alignmode | (std::max(outputformat, 1) << 8), count, refBuffer, qryBuffer, alignBuffer,
				(m_EnableBS) ? m_DirBuffer : 0);

		Log.Debug(32, "INFO\tALGN\t%d alignments computed (out of %d)", aligned, count);

		if (aligned != count)
		Log.Error("Error aligning outputs (%i of %i aligned)", aligned, count);

		//process results
		for (int i = 0; i < aligned; ++i) {
			MappedRead * cur_read = reads[i].read;
			int scoreID = reads[i].scoreId;
			int id = cur_read->ReadId;

			assert(cur_read->hasCandidates());
			cur_read->Scores[scoreID].Location.m_Location += alignBuffer[i].PositionOffset - (corridor >> 1);

			cur_read->Alignments[scoreID] = alignBuffer[i];

			Log.Debug(2048, "READ_%d\tALGN_DETAILS\tCMR_%d\t%f\t%f\t%llu\t%.*s\t%s", cur_read->ReadId, scoreID, cur_read->Scores[scoreID].Score.f, alignBuffer[i].Identity, alignBuffer[i].NM, refMaxLen, refBuffer[i], qryBuffer[i]);

			if ((cur_read->Calculated - 1) == scoreID) {

				debugAlgnFinished(cur_read);

				SaveRead(cur_read);
			}

		}
		alignTime = tmr.ET();
	} else {
		Log.Debug(1, "INFO\tALGN\tEmpty buffer submitted.");
	}
}

Align AlignmentBuffer::computeAlignment(uloc const position, int corridor,
		char * const readSeq, size_t const readLength, int const QStart,
		int const QEnd, int fullReadLength) {

	Align align;

	//initialize arrays for CIGAR and MD string
	align.pBuffer1 = new char[readLength * 4];
	align.pBuffer2 = new char[readLength * 4];
	*(int*) align.pBuffer1 = 0x212121;
	*(int*) align.pBuffer2 = 0x212121;

	int * clipping = new int[2];
	clipping[0] = QStart;
	clipping[1] = QEnd;

	int cigarLength = -1;
	//Required if corridor is increased!
	int correctOffset = (corridor >> 1);

	bool shortAlignment = false;

	int maxTries = 3;

	while (cigarLength == -1) {

		size_t const refSeqLen = readLength + corridor + 1;
		//TODO: check why decoded/SequenceProvider writes outside of refSeqLen
		//This makes + 100 necessary
		char * refSeq = new char[refSeqLen + 100];

		//decode reference sequence
		if (!SequenceProvider.DecodeRefSequenceExact(refSeq, position, refSeqLen, corridor)) {
			//Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
			Log.Warning("Could not decode reference for alignment");
			memset(refSeq, 'N', refSeqLen);
			exit(0);
		}

		//Local alignment
		if (pacbioDebug) {
			Log.Message("Aligning %d bp to %d bp (corridor %d)", readLength, refSeqLen, corridor);
			Log.Message("Additional QStart: %d, QEnd: %d", clipping[0], clipping[1]);
		}

		//test-align python (ngila)
//		SequenceLocation loc;
//		loc.m_Location = position;
//		SequenceProvider.convert(loc);
//		int len = 0;
//		printf("%llu\t%s\t%d\t%s\t%s", loc.m_Location, SequenceProvider.GetRefName(loc.getrefId(), len), corridor, refSeq, readSeq);

		if (pacbioDebug) {
			Log.Message("Ref: %s", refSeq);
			Log.Message("Read: %s", readSeq);
		}

		int mode = 0;
		cigarLength = aligner->SingleAlign(mode, corridor,
				(char const * const ) refSeq, (char const * const ) readSeq,
				align, clipping);

		int aligned = readLength - (align.QStart - clipping[0])
				- (align.QEnd - clipping[1]);

		if (aligned < 0.1f * readLength) {
			align.Score = -1.0f;
		}

		if (cigarLength == -1 /*|| aligned < 0.1f * readLength*/) {
			corridor = corridor * 2;
			correctOffset = (corridor >> 1);

			align.Score = -1.0f;
//			if (aligned < 0.1f * readLength) {
//				Log.Message("Aligned only %d out of %d bp", aligned, readLength);
//				shortAlignment = true;
//				if(maxTries-- > 0) {
//					cigarLength = -1;
//				}
//			}
			if (pacbioDebug) {
				Log.Message("Invalid alignment found. Running again with corridor %d, %d attempts left", corridor, maxTries);
				NGM.Stats->invalidAligmentCount += 1;
			}
		} else {
			if (pacbioDebug) {
				Log.Message("%d of %d bp successfully aligned", aligned, readLength);
			}
			//		getchar();
		}

		delete[] refSeq;
		refSeq = 0;
	}

	if (cigarLength != fullReadLength) {
		Log.Error("CIGAR string invalid: seq len: %d, cigar len: %d", fullReadLength,
				cigarLength);
		throw 1;
	}

//	if (shortAlignment) {
//		Log.Message("Waiting");
//		getchar();
//	}

	delete[] clipping;
	clipping = 0;

	align.PositionOffset -= correctOffset;
	return align;
}

//Align AlignmentBuffer::computeAlignment(MappedRead* read, int const scoreId,
//		int const corridor) {
//
//	Align align;
//	LocationScore & score = read->Scores[scoreId];
//
//	Log.Message("Computing alignment (%d) for position: %llu", scoreId, score.Location.m_Location);
//	Log.Message("Corridor: %d", corridor);
//	char * refBuffer = new char[read->length + corridor + 1];
//
//	//decode reference sequence
//	if (!SequenceProvider.DecodeRefSequence(refBuffer, 0,
//			score.Location.m_Location - (corridor >> 1), read->length + corridor)) {
////Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
//		Log.Warning("Could not decode reference for alignment (read: %s)", read->name);
//		memset(refBuffer, 'N', read->length * 1.2f);
//	}
//
//	//initialize arrays for CIGAR and MD string
//	align.pBuffer1 = new char[read->length * 4];
//	align.pBuffer2 = new char[read->length * 4];
//	*(int*) align.pBuffer1 = 0x212121;
//	*(int*) align.pBuffer2 = 0x212121;
//
//	int const mode = 0;
//	if (score.Location.isReverse()) {
//		aligner->SingleAlign(mode, corridor, (char const * const ) refBuffer,
//				(char const * const ) read->RevSeq, align, 0);
//		printf(">Ref_%s\n%s\n>%s_rev\n%s", read->name, refBuffer, read->name,
//				read->RevSeq);
//	} else {
//		aligner->SingleAlign(mode, corridor, (char const * const ) refBuffer,
//				(char const * const ) read->Seq, align, 0);
//		printf(">Ref_%s\n%s\n>%s\n%s", read->name, refBuffer, read->name,
//				read->Seq);
//	}
//
//	score.Location.m_Location += align.PositionOffset - (corridor >> 1);
//
//	delete[] refBuffer;
//
//	return align;
//}

#include <seqan/graph_algorithms.h>
#include <iostream>
//using namespace seqan;

using seqan::appendValue;
using seqan::String;
using seqan::Block;
using seqan::longestIncreasingSubsequence;
using seqan::length;

void printDotPlotLine(int const id, char const * const name,
		int const onReadStart, int const onReadStop, loc const onRefStart,
		loc const onRefStop, float const score, bool const isReverse,
		int const type, int const status) {
//	printf("%d\t%s\t%d\t%d\t%llu\t%llu\t%f\t%d\t%d\t%d\n", id, name,
//			onReadStart, onReadStop, onRefStart, onRefStop, score, isReverse,
//			type, status);
}

int * AlignmentBuffer::cLIS(Anchor * anchors, int const anchorsLenght,
		int & lisLength) {

	int * DP = 0;
	int * trace = 0;

	int * lis = new int[anchorsLenght];
	lisLength = 0;

	if (anchorsLenght > 0) {
		DP = new int[anchorsLenght];
		trace = new int[anchorsLenght];

		int maxLength = 1;
		int bestEnd = 0;
		DP[0] = 1;
		trace[0] = -1;

		for (int i = 0; i < anchorsLenght; i++) {
			DP[i] = 1;
			trace[i] = -1;

			for (int j = i - 1; j >= 0; j--) {
				uloc iRef =
						(anchors[i].isReverse) ?
								-1 * anchors[i].onRef : anchors[i].onRef;
				uloc jRef =
						(anchors[j].isReverse) ?
								-1 * anchors[j].onRef : anchors[j].onRef;
				uloc diff = iRef - jRef;

				if (DP[j] + 1 > DP[i] && /*jRef < iRef && */diff < 512 * 2.0f
						&& (diff > 512 * 0.5f
								|| anchors[i].onRead == anchors[j].onRead)
						&& anchors[j].isReverse == anchors[i].isReverse) {
					DP[i] = DP[j] + 1;
					trace[i] = j;
				}
			}

			if (DP[i] > maxLength) {
				bestEnd = i;
				maxLength = DP[i];
			}
		}

		while (trace[bestEnd] != -1) {
			lis[lisLength++] = bestEnd;
			bestEnd = trace[bestEnd];
		}
		lis[lisLength++] = bestEnd;

		delete[] DP;
		DP = 0;
		delete[] trace;
		trace = 0;
	}

	return lis;
}

void AlignmentBuffer::getSeqAnLIS(int allFwdAnchorsLength, int id,
		Anchor* allFwdAnchors, char* name) {

	//TODO: try longest common subsequence with 1 specific symbol per read part.
	//Second sequence consists of all matching positions on the reference (orderd)
	//Positions get the same symbol as a read part if they got a alignment score
	String<loc> seqFwd;
	String<loc> seqRev;
	String<float> weightsFwd;
	String<float> weightsRev;
	String<loc, Block<> > posFwd;
	String<loc, Block<> > posRev;
	for (int i = 0; i < allFwdAnchorsLength; ++i) {
		if (allFwdAnchors[i].isReverse) {
			appendValue(seqFwd, allFwdAnchors[i].onRef * -1);
		} else {
			appendValue(seqFwd, allFwdAnchors[i].onRef);
		}
		appendValue(weightsFwd, allFwdAnchors[i].score);
	}
	//	for (int i = 0; i < allRevAnchorsLength; ++i) {
	//		appendValue(seqRev, allRevAnchors[i].onRef);
	//		appendValue(weightsRev, allRevAnchors[i].score);
	//
	//	}
	longestIncreasingSubsequence(seqFwd, posFwd);
//	float w = heaviestIncreasingSubsequence(seqFwd, weightsFwd, posFwd);

	if (pacbioDebug)
		std::cerr << "SeqAN LIS: ";
	for (int i = length(posFwd) - 1; i >= 0; --i) {
		if (pacbioDebug)
			std::cerr << posFwd[i] << ", ";
//		printDotPlotLine(id, name, allFwdAnchors[posFwd[i]].onRead * 512,
//				allFwdAnchors[posFwd[i]].onRead * 512 + 512,
//				allFwdAnchors[posFwd[i]].onRef,
//				allFwdAnchors[posFwd[i]].onRef + 512,
//				allFwdAnchors[posFwd[i]].score,
//				allFwdAnchors[posFwd[i]].isReverse, TYPE_LIS_SEQAN, STATUS_OK);
	}
	if (pacbioDebug)
		std::cerr << std::endl;

}

bool AlignmentBuffer::isSameDirection(Interval a, Interval b) {
	return a.isReverse == b.isReverse;
}

bool AlignmentBuffer::isCompatible(Interval a, Interval b) {

//	if (!isCompatible(a,b)) {
//		return false;
//	}

	bool overlapsOnRead = false;
	bool overlapsOnRef = false;

	bool isInCorridor = false;

	double k1 = (double) (a.onReadStop - a.onReadStart)
			/ (a.onRefStop - a.onRefStart) * 1.0f;
	double d1 = (double) (a.onRefStart * 1.0f * a.onReadStop
			- a.onReadStart * 1.0f * a.onRefStop) / (a.onRefStart - a.onRefStop)
			* 1.0f;

	double k2 = (double) (b.onReadStop - b.onReadStart)
			/ (double) (b.onRefStop - b.onRefStart) * 1.0f;
	double d2 = (double) (b.onRefStart * 1.0f * b.onReadStop
			- b.onReadStart * 1.0f * b.onRefStop)
			/ (double) (b.onRefStart - b.onRefStop) * 1.0f;

//	Log.Message("a: %d, b: %d, diff: %d", d1, d2, d2 - d1);

	loc a_y0 = -(d1 / k1);
	loc b_y0 = -(d2 / k2);

	if (pacbioDebug)
		Log.Message("k1: %f, d1: %f", k1, d1);
	if (pacbioDebug)
		Log.Message("k2: %f, d2: %f", k2, d2);

	if (pacbioDebug) {
		Log.Message("Read start - a: %llu, b: %llu, diff: %llu", a_y0, b_y0, abs(b_y0 - a_y0));
	}

	isInCorridor = abs(b_y0 - a_y0) < 4096;

	return !overlapsOnRead && !overlapsOnRef && isInCorridor;
}

bool AlignmentBuffer::isContained(Interval a, Interval b) {
	return a.onReadStart >= b.onReadStart && a.onReadStop <= b.onReadStop
			&& a.onRefStart >= b.onRefStart && a.onRefStop <= b.onRefStop
			&& a.isReverse == b.isReverse;
}

bool AlignmentBuffer::isContainedOnRead(Interval a, Interval b) {
	if (abs(a.onReadStop - a.onReadStart) > abs(b.onReadStop - b.onReadStart)) {
		Interval tmp;
		tmp = a;
		a = b;
		b = tmp;
	}
	bool isContained = a.onReadStart >= b.onReadStart
			&& a.onReadStop <= b.onReadStop;

	return isContained;
}

AlignmentBuffer::Interval AlignmentBuffer::splitInterval(Interval a,
		Interval b) {

	Log.Message("=======================");
	a.print();
	b.print();
//	getchar();
	return a;
}

AlignmentBuffer::Interval AlignmentBuffer::mergeIntervals(Interval a,
		Interval b) {
	a.onReadStart = std::min(a.onReadStart, b.onReadStart);
	a.onReadStop = std::max(a.onReadStop, b.onReadStop);
	a.onRefStart = std::min(a.onRefStart, b.onRefStart);
	a.onRefStop = std::max(a.onRefStop, b.onRefStop);
	a.score += b.score;
	return a;
}

//bool AlignmentBuffer::constructMappedSegements(Interval * intervals,
//		Interval interval, int & intervalsIndex) {
//	bool addInterval = true;
//	bool done = false;
//	for (int j = 0; j < intervalsIndex && !done; ++j) {
//		if (isContained(interval, intervals[j])) {
//			if (pacbioDebug)
//				Log.Message("Interval is contained in %d", j);
//				addInterval = false;
//				done = true;
//				//Do not add
//			} else {
//				if (pacbioDebug)
//				Log.Message("Interval not contained in %d", j);
//
//				if(isCompatible(interval, intervals[j])) {
//					if (pacbioDebug)
//					Log.Message("Is compatible");
//
//					if(isSameDirection(interval, intervals[j])) {
//						addInterval = false;
//						//Merge
//						intervals[j] = mergeIntervals(interval, intervals[j]);
//						done = true;
//					} else {
//
//						if(isContainedOnRead(interval, intervals[j])) {
//							//Split interval
//							intervals[intervalsIndex++] = splitInterval(interval, intervals[j]);
//							intervals[intervalsIndex++] = interval;
//							addInterval = false;
//							done = true;
//						} else {
//							if(pacbioDebug) {
//								Log.Message("Is compatible but not contained on read - skipping");
//							}
//						}
//					}
//				} else {
//					if (pacbioDebug)
//					Log.Message("Is NOT compatible");
//					addInterval = addInterval && true;
//					//Add
//				}
//			}
//		}
//
//	if (addInterval || intervalsIndex == 0) {
//		if (intervalsIndex < 1000) {
//			intervals[intervalsIndex++] = interval;
//		} else {
//			return false;
//		}
//		if (pacbioDebug)
//			Log.Message("Adding interval");
//		}
//	return true;
//}

bool AlignmentBuffer::constructMappedSegements(MappedSegment * segments,
		Interval interval, size_t & segmentsIndex) {
//	bool addInterval = true;
//	bool done = false;
	for (int i = 0; i < segmentsIndex; ++i) {
		Interval * intervals = segments[i].list;
		for (int j = 0; j < segments[i].length; ++j) {
			if (isContained(interval, intervals[j])) {
				if (pacbioDebug) {
					Log.Message("Interval is contained in %d", j);
				}

				//Do not add
				return true;
			} else {
				if (pacbioDebug) {
					Log.Message("Interval not contained in %d", j);
				}

				if(isCompatible(interval, intervals[j])) {
					if (pacbioDebug) {
						Log.Message("Adding interval to segment");
					}
					if (segments[i].length < 100) {
						intervals[segments[i].length++] = interval;
						return true;
					} else {
						return false;
					}
				} else {
					if (pacbioDebug) {
						Log.Message("Is NOT compatible");
					}
					//Add
				}
			}
		}
	}

	if (pacbioDebug) {
		Log.Message("Creating new segment from interval");
	}
	segments[segmentsIndex].list[0] = interval;
	segments[segmentsIndex++].length = 1;

	return true;
}

void AlignmentBuffer::consolidateSegment(Interval * intervals,
		int & intervalsIndex, MappedSegment segment) {

	Interval lastInterval = segment.list[0];
	Interval currentInterval;

	for (int i = 1; i < segment.length; ++i) {
		currentInterval = segment.list[i];
		if (isSameDirection(currentInterval, lastInterval)) {
			lastInterval = mergeIntervals(currentInterval, lastInterval);
		} else {
			//Add lastInterval
			intervals[intervalsIndex++] = lastInterval;
			lastInterval = currentInterval;
		}
	}
	intervals[intervalsIndex++] = lastInterval;

}

bool sortIntervalsInSegment(AlignmentBuffer::Interval a,
		AlignmentBuffer::Interval b) {
	return a.onReadStart < b.onReadStart;
}

AlignmentBuffer::Interval * AlignmentBuffer::consolidateSegments(
		MappedSegment * segments, size_t segmentsIndex, int & intervalsIndex) {

	if (pacbioDebug) {
		Log.Message("============================");
	}
	Interval * intervals = new Interval[1000];

	for (int i = 0; i < segmentsIndex; ++i) {

		std::sort(segments[i].list, segments[i].list + segments[i].length,
				sortIntervalsInSegment);

		if (pacbioDebug) {
			Log.Message("Segment %d: ", i);
			for(int j = 0; j < segments[i].length; ++j) {
				Log.Message("\tInterval %d: ", j);
				segments[i].list[j].print();
			}
		}

		consolidateSegment(intervals, intervalsIndex, segments[i]);
	}
	if (pacbioDebug) {
		Log.Message("============================");
	}

	return intervals;
}

AlignmentBuffer::Interval * AlignmentBuffer::findSubsequences(char * name,
		int id, Anchor * allFwdAnchors, int allFwdAnchorsLength,
		Anchor * allRevAnchors, int allRevAnchorsLength, int readParts,
		int readLenth, int & intervalsIndex) {

	if (pacbioDebug) {
		Log.Message("Finding LIS for read %d", id);
	}

	//Interval * intervals = new Interval[1000];
	MappedSegment * segments = new MappedSegment[100];
	size_t segementsIndex = 0;
	//intervalsIndex = 0;

	//Sort by position on read. Prpbably not necessary!!
	std::sort(allFwdAnchors, allFwdAnchors + allFwdAnchorsLength,
			AlignmentBuffer::sortAnchorOnRead);

	//Get standard LIS with Seqan
	getSeqAnLIS(allFwdAnchorsLength, id,
			allFwdAnchors, name);

	int run = 0;
	bool finished = false;
	while (!finished) {
		if (allFwdAnchorsLength > 0) {

			//Print remaining elements
			if (pacbioDebug) {
				std::cerr << "Elements: ";
				for (int i = 0; i < allFwdAnchorsLength; ++i) {
					std::cerr << i << ":" << allFwdAnchors[i].onRead << ":"
					<< allFwdAnchors[i].onRef << ":"
					<< allFwdAnchors[i].isReverse << ", ";
				}
				std::cerr << std::endl;
			}

			//Find constrained LIS
			int lisLength = 0;
			int * lis = cLIS(allFwdAnchors, allFwdAnchorsLength, lisLength);

			int minOnRead = 999999;
			int maxOnRead = 0;

			loc minOnRef = 999999999999;
			loc maxOnRef = 0;

			bool isReverse = false;

			float intervalScore = 0.0f;

			if (pacbioDebug)
			std::cerr << "cLIS" << run << ":     ";
			//Remove LIS from candidates
			int posInLIS = lisLength - 1;
			allRevAnchorsLength = 0;
			for (int i = 0; i < allFwdAnchorsLength; ++i) {
				if (posInLIS >= 0 && i == lis[posInLIS]) {

					int onRead = allFwdAnchors[lis[posInLIS]].onRead * 512;
//					if(allFwdAnchors[lis[posInLIS]].isReverse) {
//						onRead = readLenth - onRead;
//					}

					//Print current LIS
					if (pacbioDebug)
					std::cerr << lis[posInLIS] << ", ";
					printDotPlotLine(id, name,
							onRead,
							onRead + 512,
							allFwdAnchors[lis[posInLIS]].onRef,
							allFwdAnchors[lis[posInLIS]].onRef + 512,
							allFwdAnchors[lis[posInLIS]].score,
							allFwdAnchors[lis[posInLIS]].isReverse,
							TYPE_LIS_SEQAN + run, STATUS_OK);

					minOnRead = std::min(minOnRead, onRead);
					maxOnRead = std::max(maxOnRead, onRead + 512);

					minOnRef = std::min(minOnRef, allFwdAnchors[lis[posInLIS]].onRef);
					maxOnRef = std::max(maxOnRef, allFwdAnchors[lis[posInLIS]].onRef + 512);

					isReverse = allFwdAnchors[lis[posInLIS]].isReverse;
					intervalScore += allFwdAnchors[lis[posInLIS]].score;

					//Remove from remaining elements
					posInLIS -= 1;
				} else {
					//Add to second array
					allRevAnchors[allRevAnchorsLength++] = allFwdAnchors[i];
					//			printDotPlotLine(id, name, allFwdAnchors[i].onRead,
					//					allFwdAnchors[i].onRef, allFwdAnchors[i].score,
					//					allFwdAnchors[i].isReverse, TYPE_CLIS_2, STATUS_OK);
				}
			}

			if (pacbioDebug) {
				std::cerr << std::endl;
			}

			Interval interval;

			int addStart = std::min(512, minOnRead);
			int addEnd = std::min(512, readLenth - maxOnRead);

			interval.onReadStart = minOnRead - addStart;
			interval.onReadStop = maxOnRead + addEnd;

			interval.onRefStart = minOnRef - addStart;
			interval.onRefStop = maxOnRef + addEnd;

			interval.isReverse = isReverse;
			interval.score = intervalScore;

			printDotPlotLine(id, name,
					interval.onReadStart,
					interval.onReadStop,
					interval.onRefStart,
					interval.onRefStop,
					interval.score,
					interval.isReverse,
					TYPE_CLIS + run, STATUS_OK);

			//TODO: check chromosome borders

			if (pacbioDebug) {
				interval.print();
			}

			if(!constructMappedSegements(segments, interval, segementsIndex)) {
				Log.Message("Too many intervals found (>1000). Ignoring all others");
				lisLength = 0;
				delete[] lis;
//					return intervals;
				Fatal();
			}
			if(pacbioDebug) {
				Log.Message("%d segments in list", segementsIndex);
			}

			//Switch first with second list
			Anchor * tmp = allFwdAnchors;
			allFwdAnchors = allRevAnchors;
			allFwdAnchorsLength = allRevAnchorsLength;
			allRevAnchors = tmp;
			allRevAnchorsLength = 0;

			lisLength = 0;
			delete[] lis;

			run += 1;
			if (run == 4) {
				finished = true;
			}
		} else {
			finished = true;
		}
	}

	intervalsIndex = 0;
	Interval * intervals = consolidateSegments(segments, segementsIndex, intervalsIndex);
	//delete intervals;
	//intervals = 0;
	delete[] segments;
	segments = 0;

	return intervals;
}

static inline char cplBase(char c) {
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

void computeReverseSeq(char * Seq, char * RevSeq, int const length) {
	//static int const qryMaxLen = Config.GetInt("qry_max_len");
//	RevSeq = new char[length + 1];
//	memset(RevSeq, 0, length + 1);

	char * fwd = Seq;
	char * rev = RevSeq + length - 1;

	for (int i = 0; i < length; ++i) {
		*rev-- = cplBase(*fwd++);
	}
}

bool AlignmentBuffer::inversionDetection(Align const align,
		Interval const interval, int const readLength, char * fullReadSeq,
		Interval & leftOfInv, Interval & rightOfInv, Interval & inv,
		char const * const readName) {

	static bool const enabled = Config.GetInt(NOINVERSIONS) != 1;

	if (!enabled) {
		return false;
	}

	loc start = interval.onRefStart;

	SequenceLocation seqLoc;
	seqLoc.m_Location = start + align.PositionOffset;

	if (pacbioDebug) {
		SequenceProvider.convert(seqLoc);
		Log.Message("Mapping Pos: %llu", seqLoc.m_Location);
	}

	bool inversionFound = false;

	int * extData = (int *) align.ExtendedData;
	int i = 0;
	int len = 0;

	int alignStartOnRef = extData[i++];
	int alignStartOnRead = extData[i++];

	int alignStopOnRef = extData[i++];
	int alignStopOnRead = extData[i++];

	while (extData[i] > 0 && !inversionFound) {
		uloc ref = extData[i++];
		uloc read = extData[i++];

		char * refSeq = new char[250];
		char * readSeq = new char[500];
		char * revreadSeq = new char[500];
		memset(refSeq, '\0', 250 * sizeof(char));
		memset(readSeq, '\0', 500 * sizeof(char));
		memset(revreadSeq, '\0', 500 * sizeof(char));

		SequenceLocation testLoc;
		testLoc.m_Location = start + align.PositionOffset + ref - 100;
		SequenceProvider.DecodeRefSequence(refSeq, 0, testLoc.m_Location, 200);
		strncpy(readSeq, fullReadSeq + read - 200, 400);

		computeReverseSeq(readSeq, revreadSeq, 400);

		if (pacbioDebug) {
			Log.Message("Ref:     %s", refSeq);
			Log.Message("Read:    %s", readSeq);
			Log.Message("RevRead: %s", revreadSeq);
		}

		IAlignment * aligner = new EndToEndAffine();

		float scoreFwd = 0.0f;
		float scoreRev = 0.0f;
		aligner->SingleScore(10, 0, readSeq, refSeq, scoreFwd, 0);
		aligner->SingleScore(10, 0, revreadSeq, refSeq, scoreRev, 0);

		if (pacbioDebug) {
			Log.Message("Fwd: %f, Rev: %f", scoreFwd, scoreRev);
			SequenceProvider.convert(testLoc);
		}

//		printf("%s\t%llu\t%llu\t%s\t%d\n", SequenceProvider.GetRefName(testLoc.getrefId(), len), testLoc.m_Location, testLoc.m_Location + 200, readName, 100);
		if (scoreRev > scoreFwd) {

			inversionFound = true;
			if (pacbioDebug) {
				Log.Message("Inversion detected: %llu, %llu", ref,
						read);
				Log.Message("READ LENGTH: %d", readLength);
//				printf("%s\t%llu\t%llu\t%s\t%d\n", SequenceProvider.GetRefName(testLoc.getrefId(), len), testLoc.m_Location, testLoc.m_Location + 200, readName, 100);
			}

			if (interval.isReverse) {

				leftOfInv.onReadStop = readLength
						- (interval.onReadStart + alignStartOnRead);

				leftOfInv.onRefStart = start + align.PositionOffset
						+ alignStartOnRef;

				leftOfInv.onReadStart = readLength
						- (interval.onReadStart + read);
				leftOfInv.onRefStop = start + align.PositionOffset + ref;

				leftOfInv.isReverse = interval.isReverse;

				rightOfInv.onReadStop = readLength
						- (interval.onReadStart + read);
				rightOfInv.onRefStart = start + align.PositionOffset + ref;
				rightOfInv.onReadStart = readLength
						- (interval.onReadStart + alignStopOnRead);
				rightOfInv.onRefStop = start + align.PositionOffset
						+ alignStopOnRef;

				rightOfInv.isReverse = interval.isReverse;

				inv.onReadStart = leftOfInv.onReadStart;
				inv.onReadStop = leftOfInv.onReadStart;
				inv.onRefStart = leftOfInv.onRefStop;
				inv.onRefStop = leftOfInv.onRefStop;
				inv.isReverse = !rightOfInv.isReverse;
			} else {

				leftOfInv.onReadStart = interval.onReadStart + alignStartOnRead;
				leftOfInv.onRefStart = start + align.PositionOffset
						+ alignStartOnRef;
				leftOfInv.isReverse = interval.isReverse;
				leftOfInv.onReadStop = interval.onReadStart + read;
				leftOfInv.onRefStop = start + align.PositionOffset + ref;

				rightOfInv.onReadStart = interval.onReadStart + read;
				rightOfInv.onRefStart = start + align.PositionOffset + ref;
				rightOfInv.isReverse = interval.isReverse;
				rightOfInv.onReadStop = interval.onReadStart + alignStopOnRead;
				rightOfInv.onRefStop = start + align.PositionOffset
						+ alignStopOnRef;

				inv.onReadStart = leftOfInv.onReadStop;
				inv.onReadStop = leftOfInv.onReadStop;
				inv.onRefStart = leftOfInv.onRefStop;
				inv.onRefStop = leftOfInv.onRefStop;
				inv.isReverse = !rightOfInv.isReverse;
			}

			if (pacbioDebug) {
				SequenceLocation seqLoc2;
				seqLoc2.m_Location = leftOfInv.onRefStart;
				SequenceProvider.convert(seqLoc2);
				Log.Message("Mapping Pos2: %llu - %llu (alignStartOnREf: %llu, corridor: %d)", leftOfInv.onRefStart, seqLoc2.m_Location, alignStartOnRef, corridor);
				SequenceLocation seqLoc5;
				seqLoc5.m_Location = rightOfInv.onRefStart;
				SequenceProvider.convert(seqLoc5);
				Log.Message("Mapping Pos5: %llu - %llu", rightOfInv.onRefStart, seqLoc5.m_Location);
				SequenceLocation seqLoc3;
				seqLoc3.m_Location = rightOfInv.onRefStart;
				SequenceProvider.convert(seqLoc3);
				Log.Message("Mapping Pos3: %llu - %llu", rightOfInv.onRefStart, seqLoc3.m_Location);
				SequenceLocation seqLoc4;
				seqLoc4.m_Location = rightOfInv.onRefStop;
				SequenceProvider.convert(seqLoc4);
				Log.Message("Mapping Pos4: %llu - %llu", rightOfInv.onRefStop, seqLoc4.m_Location);
			}

		} else {
//			printf("%s\t%llu\t%llu\t%s\t%d\n", SequenceProvider.GetRefName(testLoc.getrefId(), len), testLoc.m_Location, testLoc.m_Location + 200, readName, 0);
		}

		delete[] refSeq;
		refSeq = 0;
		delete[] readSeq;
		readSeq = 0;
		delete[] revreadSeq;
		revreadSeq = 0;
		delete aligner;
		aligner = 0;
	}

	return inversionFound;
}

int AlignmentBuffer::estimateCorridor(const Interval & interval) {
	int corridor = 8192;
	int onRead = interval.onReadStop - interval.onReadStart;
	int onRef = interval.onRefStop - interval.onRefStart;
	int diff = onRead - onRef;
	int corridorFromDiff = (int) ((abs(diff) * 2.1f));
	int corridorFromLength = (int) ((abs(onRead) * 0.20f));
	corridor = std::min(corridor,
			std::max(corridorFromDiff, corridorFromLength));
	if (pacbioDebug)
		Log.Message("Corridor estimation - onRead: %d, onRef: %d, diff: %d, corridor: %d", onRead, onRef, diff, corridor);

	return corridor;
}

Align AlignmentBuffer::alignInterval(MappedRead const * const read_,
		Interval const interval, char * const readSeq,
		size_t const readSeqLen) {

	//int corridor = std::max(abs(difference) * 2, (int) (read->length * 0.2));

	Align align;
	align.Score = 0.0f;

	int corridor = estimateCorridor(interval);

	int QStart = 0;
	int QEnd = 0;

	if (interval.isReverse) {
		QEnd = interval.onReadStart;
		QStart = read_->length - interval.onReadStop;

	} else {
		QStart = interval.onReadStart;
		QEnd = read_->length - interval.onReadStop;
	}

	if (pacbioDebug) {
		Log.Message("Start pos: %d, Length: %d", interval.onReadStart, readSeqLen);
	}

	//test-align python (ngila)
//	printf("%s\t%d\t", read_->name, 0);
	if (pacbioDebug) {
		Log.Message("Computing alignment");
	}

	try {
		align = computeAlignment(interval.onRefStart, corridor, readSeq,
				readSeqLen, QStart, QEnd, read_->length);
	} catch (int e) {
		Log.Error("Error occurred while aligning read %s", read_->name);
		Fatal();
	}

	return align;
}

char* const AlignmentBuffer::extractReadSeq(const size_t& readSeqLen,
		Interval& interval, MappedRead* read) {
	if (pacbioDebug) {
		Log.Message("Allocating %d bytes", readSeqLen + 1);
	}

	char*const readSeq = new char[readSeqLen + 1];
	if (interval.isReverse) {
		read->computeReverseSeq();
		computeReverseSeq(read->Seq + interval.onReadStart, readSeq, readSeqLen);
		readSeq[readSeqLen] = '\0';
	} else {
		strncpy(readSeq, read->Seq + interval.onReadStart, readSeqLen);
		readSeq[readSeqLen] = '\0';
	}
	if (pacbioDebug) {
		Log.Message("ReadSeqLen: %d", readSeqLen);
	}
	if (pacbioDebug) {
		Log.Message("Start pos: %d, Length: %d", interval.onReadStart, readSeqLen);
	}
	return readSeq;
}

bool AlignmentBuffer::alignInversion(Interval interval, Interval leftOfInv,
		Interval inv, Interval rightOfInv, MappedRead * read, Align * tmpAling,
		int & alignIndex, LocationScore * tmp) {
	int inversionLength = 0;
	if (pacbioDebug) {
		Log.Message("Inversion detected!");
		Log.Message("Interval:");
		interval.print();
		Log.Message("Left Interval:");
		leftOfInv.print();
	}
	/**********************************************/
	//Align part of read that is left of inversion
	/**********************************************/
	size_t readSeqLen = leftOfInv.onReadStop - leftOfInv.onReadStart;
	char* readSeq = extractReadSeq(readSeqLen, leftOfInv, read);
	Align align = alignInterval(read, leftOfInv, readSeq, readSeqLen);
	/**********************************************/
	//QEnd from first read part, gives length of
	//inversion that was covered by left part
	/**********************************************/
	inversionLength += (readSeqLen - ((int*) (align.ExtendedData))[3]);
	if (inv.isReverse) {
		inv.onReadStart -= (inversionLength * 1.2);
	} else {
		inv.onReadStop += (inversionLength * 1.2);
	}
	inv.onRefStart -= (inversionLength * 1.5);
	if (pacbioDebug) {
		Log.Message("1. pos in read end: %d", inversionLength);
	}
	delete[] readSeq;
	readSeq = 0;
	delete[] (int*) (align.ExtendedData);
	align.ExtendedData = 0;
	if (align.Score > 0.0f) {
		tmpAling[alignIndex] = align;

		tmp[alignIndex].Location.m_Location = leftOfInv.onRefStart
				+ align.PositionOffset;	//- (corridor >> 1); handled in computeAlingment
		tmp[alignIndex].Location.setReverse(leftOfInv.isReverse);
		tmp[alignIndex].Score.f = align.Score;

		read->Calculated += 1;
		alignIndex += 1;

	} else {
		if (pacbioDebug) {
			Log.Message("Alignment failed");
		}
		return false;
	}
	if (pacbioDebug) {
		Log.Message("Right Interval:");
		rightOfInv.print();
	}
	/**********************************************/
	//Align part of read that is right of inversion
	/**********************************************/
	readSeqLen = rightOfInv.onReadStop - rightOfInv.onReadStart;
	readSeq = extractReadSeq(readSeqLen, rightOfInv, read);
	align = alignInterval(read, rightOfInv, readSeq, readSeqLen);
	/**********************************************/
	//QStart from second read part, gives length
	//of inversion that was covered by second read
	/**********************************************/
	inversionLength += ((int*) (align.ExtendedData))[1];
	if (inv.isReverse) {
		inv.onReadStop += (((int*) (align.ExtendedData))[1] * 1.2);
	} else {
		inv.onReadStart -= (((int*) (align.ExtendedData))[1] * 1.2);
	}
	inv.onRefStop += (((int*) (align.ExtendedData))[1] * 1.5);
	if (pacbioDebug) {
		Log.Message("2. pos in read start: %d", ((int *) align.ExtendedData)[1]);
	}
	delete[] readSeq;
	readSeq = 0;
	delete[] (int*) (align.ExtendedData);
	align.ExtendedData = 0;
	if (align.Score > 0.0f) {
		tmpAling[alignIndex] = align;

		tmp[alignIndex].Location.m_Location = rightOfInv.onRefStart
				+ align.PositionOffset;	//- (corridor >> 1); handled in computeAlingment
		tmp[alignIndex].Location.setReverse(rightOfInv.isReverse);
		tmp[alignIndex].Score.f = align.Score;

		read->Calculated += 1;
		alignIndex += 1;

	} else {
		if (pacbioDebug) {
			Log.Message("Alignment failed");
		}
		return false;
	}
	/**********************************************/
	//Align inverted part, length was inferred
	//from left and right part
	/**********************************************/
	if (pacbioDebug) {
		Log.Message("Lenght of inversion is %d", inversionLength);
		Log.Message("Inverted Interval:");
		inv.print();
	}
	//Align part of read that is right of inversion
	readSeqLen = inv.onReadStop - inv.onReadStart;
	readSeq = extractReadSeq(readSeqLen, inv, read);
	align = alignInterval(read, inv, readSeq, readSeqLen);
	delete[] readSeq;
	readSeq = 0;
	delete[] (int*) (align.ExtendedData);
	align.ExtendedData = 0;
	if (align.Score > 0.0f) {
		tmpAling[alignIndex] = align;

		tmp[alignIndex].Location.m_Location = inv.onRefStart
				+ align.PositionOffset;	//- (corridor >> 1); handled in computeAlingment
		tmp[alignIndex].Location.setReverse(inv.isReverse);
		tmp[alignIndex].Score.f = align.Score;

		read->Calculated += 1;
		alignIndex += 1;

	} else {
		if (pacbioDebug) {
			Log.Message("Alignment failed");
		}
	}

	return true;
}

void AlignmentBuffer::alignSingleOrMultipleIntervals(MappedRead * read,
		Interval interval, LocationScore * tmp, Align * tmpAling,
		int & alignIndex) {

	size_t const readSeqLen = interval.onReadStop - interval.onReadStart;
	char * readSeq = extractReadSeq(readSeqLen, interval, read);

	Align align = alignInterval(read, interval, readSeq, readSeqLen);

	if (pacbioDebug) {
		SequenceLocation seqLoc2;
		seqLoc2.m_Location = interval.onRefStart + align.PositionOffset;
		SequenceProvider.convert(seqLoc2);
		Log.Message("Mapping Pos (%s): %llu -- %llu", read->name, interval.onRefStart + align.PositionOffset, seqLoc2.m_Location);
	}

	Interval leftOfInv;
	Interval rightOfInv;
	Interval inv;
	bool inversionDetected = false;
	bool inversionAligned = false;
	if ((inversionDetected = inversionDetection(align, interval,
			read->length /*+ readSeqLen*/, readSeq, leftOfInv, rightOfInv, inv,
			read->name))) {

		inversionAligned = alignInversion(interval, leftOfInv, inv, rightOfInv,
				read, tmpAling, alignIndex, tmp);
		if (!inversionAligned) {
			SequenceLocation loc;
			loc.m_Location = inv.onRefStart;
			SequenceProvider.convert(loc);
			SequenceLocation loc2;
			loc2.m_Location = inv.onRefStop;
			SequenceProvider.convert(loc2);
			int len = 0;
			Log.Message("Potential inversion detected at %s:%llu-%llu but not taken into account due to errors. (read: %s)", SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location, loc2.m_Location, read->name);
		}
	} else {
		Log.Message("No inversion detected!");
	}

	if (!inversionDetected || !inversionAligned) {
		/**********************************************/
		// No inversion detected
		/**********************************************/
		if (align.Score > 0.0f) {
			tmpAling[alignIndex] = align;

			tmp[alignIndex].Location.m_Location = interval.onRefStart
					+ align.PositionOffset;	//- (corridor >> 1); handled in computeAlingment
			tmp[alignIndex].Location.setReverse(interval.isReverse);
			tmp[alignIndex].Score.f = align.Score;

			read->Calculated += 1;
			alignIndex += 1;
		} else {
			if (pacbioDebug) {
				Log.Message("Alignment failed");
			}

		}

		/**********************************************/
		// Check if significant portion of read was
		// not aligned
		/**********************************************/
		int QStart = ((int*) (align.ExtendedData))[1];
		int QEnd = readSeqLen - ((int*) (align.ExtendedData))[3];

//	if (QStart > 0.1f * readSeqLen && QStart > 512) {
//
//		SequenceLocation loc;
//		loc.m_Location = tmp[alignIndex - 1].Location.m_Location;
//		SequenceProvider.convert(loc);
//		int len = 0;
//		Log.Message("Read: %s (%s:%llu)", read->name, SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location);
//		Log.Message("Skipped %d bases at beginning of %d bp long read (%f)", QStart, readSeqLen, (QStart * 100.0f / readSeqLen));
//		Log.Message("%s", align.pBuffer1);
//		Log.Message("Waiting");
//		getchar();
//	}
//
//	if (QEnd > 0.1f * readSeqLen && QStart > 512) {
//		SequenceLocation loc;
//		loc.m_Location = tmp[alignIndex - 1].Location.m_Location;
//		SequenceProvider.convert(loc);
//		int len = 0;
//		Log.Message("Read: %s (%s:%llu)", read->name, SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location);
//		Log.Message("Skipped %d bases at end of %d bp long read (%f)", QEnd, readSeqLen, (QEnd * 100.0f / readSeqLen));
//		Log.Message("%s", align.pBuffer1);
//		Log.Message("Waiting");
//		getchar();
//	}
	}
	delete[] (int *) align.ExtendedData;
	align.ExtendedData = 0;
	delete[] readSeq;
	readSeq = 0;

}

void reconcileRead(ReadGroup * group) {

	for (int j = 0; j < group->readNumber; ++j) {
		MappedRead * part = group->reads[j];

		Log.Message("%d: %d", j, part->mappingQlty);
	}

	MappedRead * read = group->fullRead;

	for (int i = 0; i < read->Calculated; ++i) {
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%f\n", read->name, i, read->length,
				read->Alignments[i].QStart, read->Alignments[i].QEnd,
				read->Scores[i].Location.isReverse(),
				read->Alignments[i].Score);
	}
	getchar();
}

void AlignmentBuffer::processLongReadLIS(ReadGroup * group) {

	Timer tmr;
	tmr.ST();

	if (pacbioDebug) {
		Log.Message("Processing LongReadLIS: %d - %s (lenght %d)", group->fullRead->ReadId, group->fullRead->name, group->fullRead->length);
		Log.Message("Read length: %d", group->fullRead->length);
	}

	float avgGroupScore = group->bestScoreSum * 1.0f / group->readsFinished;
	float minGroupScore = avgGroupScore * 0.5f;

//	bool isReverse = false;
	if (group->fwdMapped < group->reverseMapped) {
		//Read most probably maps on reverse strand
		if (pacbioDebug) {
			Log.Message("Read is reverse!");
		}
//		isReverse = true;
	}

	Anchor * anchorsFwd = new Anchor[100000];
	int anchorFwdIndex = 0;

	Anchor * anchorsRev = new Anchor[100000];
	int anchorRevIndex = 0;

	for (int j = 0; j < group->readNumber; ++j) {
		MappedRead * part = group->reads[j];

		//Reverse order if read is reverse
//		int onRead = (isReverse) ? (int) (group->readNumber) - j : j;
		int onRead = j;

		float minScore =
				(part->numScores() > 0) ? part->Scores[0].Score.f * 0.8 : 0.0f;
		minScore = std::max(minScore, minGroupScore);
//		minScore = 0.0f;
		bool print = false;

		int maxCandidatesPerReadPart = 100;
		for (int k = 0; k < part->numScores(); ++k) {
			//Since order of read parts is reversed for reverse reads, reverse flag for the parts has to be changed as well
//			bool scoreReverse =
//			(isReverse) ?
//			!part->Scores[k].Location.isReverse() :
//			part->Scores[k].Location.isReverse();
			bool scoreReverse = part->Scores[k].Location.isReverse();
			if (part->numScores() < maxCandidatesPerReadPart) {
				if (part->Scores[k].Score.f > minScore) {
//					Anchor & anchor =
//							(part->Scores[k].Location.isReverse()) ?
//									anchorsRev[anchorRevIndex++] :
//									anchorsFwd[anchorFwdIndex++];
					Anchor & anchor = anchorsFwd[anchorFwdIndex++];

					anchor.onRead = onRead;
					anchor.onRef = part->Scores[k].Location.m_Location;
					anchor.score = part->Scores[k].Score.f;
					anchor.isReverse = scoreReverse;
					anchor.type = STATUS_OK;
					//Log.Message("%s\t%d\t%llu\t%f", group->fullRead->name, j, part->Scores[k].Location.m_Location, part->Scores[k].Score.f);
					if (pacbioDebug)
						Log.Message("\t%d\t%f at %llu", j, part->Scores[k].Score.f, part->Scores[k].Location.m_Location);

					printDotPlotLine(group->fullRead->ReadId,
							group->fullRead->name, anchor.onRead * 512,
							anchor.onRead * 512 + 512,
							part->Scores[k].Location.m_Location,
							part->Scores[k].Location.m_Location + 512,
							part->Scores[k].Score.f,
							part->Scores[k].Location.isReverse(),
							TYPE_UNFILTERED, STATUS_OK);
					print = true;

				} else {

					//Score too low
					printDotPlotLine(group->fullRead->ReadId,
							group->fullRead->name, onRead * 512,
							onRead * 512 + 512,
							part->Scores[k].Location.m_Location,
							part->Scores[k].Location.m_Location + 512,
							part->Scores[k].Score.f,
							part->Scores[k].Location.isReverse(),
							TYPE_UNFILTERED, STATUS_LOWSCORE);
					print = true;
				}
			} else {
				//Repetitive
				printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
						onRead * 512, onRead * 512 + 512,
						part->Scores[k].Location.m_Location,
						part->Scores[k].Location.m_Location + 512,
						part->Scores[k].Score.f,
						part->Scores[k].Location.isReverse(), TYPE_UNFILTERED,
						STATUS_REPETITIVE);
				print = true;
			}
		}

		if (!print) {
			//No hits found
			if (pacbioDebug) {
				Log.Message("No hits found for part %d", onRead);
				printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
						onRead * 512, onRead * 512 + 512, 0, 0, 0.0f, 0, TYPE_UNFILTERED, STATUS_NOHIT);
			}
		}
	}

	int intervalsIndex = 0;
	Interval * intervals = findSubsequences(group->fullRead->name,
			group->fullRead->ReadId, anchorsFwd, anchorFwdIndex, anchorsRev,
			anchorRevIndex, group->readNumber, group->fullRead->length,
			intervalsIndex);

	if (intervalsIndex != 0) {

		if (pacbioDebug) {
			Log.Message("================Intervalls found================");
			for (int i = 0; i < intervalsIndex; ++i) {
				Interval interval = intervals[i];
				interval.printOneLine();
			}
			Log.Message("================++++++++++++++++================");
		}

		//Prepare alignment list in read object
		MappedRead * read = group->fullRead;

		LocationScore * tmp = new LocationScore[intervalsIndex * 4];
		Align * tmpAling = new Align[intervalsIndex * 4];
		int alignIndex = 0;

		read->Calculated = 0;

		if (pacbioDebug)
		Log.Message("========================");
		Timer tmr;
		if (pacbioDebug) {
			tmr.ST();
		}
		for (int i = 0; i < intervalsIndex; ++i) {

			if (pacbioDebug) {
				intervals[i].print();
				intervals[i].print(group->fullRead->ReadId, group->fullRead->name,
						i);
			}

			alignSingleOrMultipleIntervals(read, intervals[i], tmp, tmpAling, alignIndex);

		}
		if (pacbioDebug) {
			Log.Message("Alignment took %fs", tmr.ET());
		}

		if (pacbioDebug) {
			Log.Message("================Intervalls aligned================");
			int bpAligned = 0;

			for(int i = 0; i < read->Calculated; ++i) {
				Align align = tmpAling[i];
				Log.Message("%d - %d", align.QStart, read->length - align.QEnd);
				bpAligned += (read->length - align.QStart - align.QEnd);
			}
			Log.Message("Aligned %d bp of %d bp", bpAligned, read->length);
			Log.Message("================++++++++++++++++++================");
		}

		read->AllocScores(tmp, alignIndex);
		read->Alignments = tmpAling;

		delete[] tmp;
		tmp = 0;

		if (pacbioDebug) {
			Log.Message("========================");
		}
		if (read->Calculated > 0) {

			reconcileRead(group);

			WriteRead(group->fullRead, true);
		} else {
			WriteRead(group->fullRead, false);
		}
	} else {
		//No candidates found for read
		WriteRead(group->fullRead, false);
	}
//	getchar();

	delete[] anchorsFwd;
	anchorsFwd = 0;
	delete[] anchorsRev;
	anchorsRev = 0;
	delete[] intervals;
	intervals = 0;

	alignTime += tmr.ET();
}

void AlignmentBuffer::processLongRead(ReadGroup * group) {
//					Log.Message("Read group with id %d finished", group->readId);
//					Log.Message("Name: %s", cur_read->name);
//					Log.Message("Reads in group: %d", group->readNumber);
//					Log.Message("Reads finished: %d", group->readsFinished);
//					Log.Message("Fwd: %d, Rev: %d", group->fwdMapped, group->reverseMapped);
//					Log.Message("Avg best score %f", group->bestScoreSum * 1.0f / group->readsFinished);

	Log.Error("Should not be here. AlignmentBuffer processLongRead");
	Fatal();

	float avgGroupScore = group->bestScoreSum * 1.0f / group->readsFinished;
	float minGroupScore = avgGroupScore * 0.8f;

//	for (int j = 0; j < group->readNumber; ++j) {
//		MappedRead * part = group->reads[j];
//		//Log.Message("ID: %d (has %d scores)", part->ReadId, part->numScores());
//		float minScore = part->Scores[0].Score.f * 0.8;
//		for (int k = 0; k < part->numScores(); ++k) {
//			if (part->Scores[k].Score.f > minScore) {
//				//Log.Message("\t%f at %llu", part->Scores[k].Score.f, part->Scores[k].Location.m_Location);
//			}
//		}
//	}

//Find first read part that maps with a min score
	int first = 0;
	while (first < group->readNumber
			&& (group->reads[first]->numScores() == 0
					|| group->reads[first]->Scores[0].Score.f < minGroupScore)) {
		first += 1;
	}

//Find last read part that maps with a min score
	int last = group->readNumber - 1;
	while (last >= 0
			&& (group->reads[last]->numScores() == 0
					|| group->reads[last]->Scores[0].Score.f < minGroupScore)) {
		last -= 1;
	}

	if (first == group->readNumber || last < 0 || group->readNumber < 3) {
		Log.Verbose("Could not map read.");
	} else {
//Distance on read between start of first and last mapped read part
//+1 to take the full last read part into account
		int distOnRead = (last - first + 1) * 512;

//If not the whole read is aligned add half of the read part size to alignment
		size_t startPosOnRead = std::max(0, first * 512 - 256);
		size_t endPosOnRead = std::min(group->fullRead->length, (last + 1) * 512 + 256);

		if (pacbioDebug)
		Log.Message("Name: %s", group->fullRead->name);
		if (pacbioDebug)
		Log.Message("On read (length %d): %d to %d (dist %d)", group->fullRead->length, startPosOnRead, endPosOnRead, distOnRead);

		bool isReverse = false;
//Compute distance between mapped location of first and last mapped read part on reference
		uloc endPos = group->reads[last]->Scores[0].Location.m_Location;
		uloc startPos = group->reads[first]->Scores[0].Location.m_Location;
		int distOnRef = 0;
		if(startPos > endPos) {
			distOnRef = startPos - endPos;
			isReverse = true;
		} else {
			distOnRef = endPos - startPos;
		}
//+512 to take full length of last read part into account
		distOnRef += 512;
		if (pacbioDebug)
		Log.Message("Start pos on ref: %llu to %llu (dist %llu)", startPos, endPos, distOnRef);

//Log.Message("Start: %llu, End: %llu", first, last);
//Log.Message("Read length: %d", group->fullRead->length);

		float coveredOnRead = (distOnRead) * 100.0f / group->fullRead->length;
		if (pacbioDebug)
		Log.Message("Covered on read: %f", coveredOnRead);
//Log.Message("On read: %d, On ref: %llu", distOnRead, distOnRef);

//Difference between distance in read cooridnates and distance in ref coordinates
//If read doesn't span larger structural variations, difference should only be
//caused by PacBio sequence error model
		int difference = distOnRead - distOnRef;
		float diffPerc = (distOnRead - distOnRef) * 1.0f / group->fullRead->length;
		if (pacbioDebug)
		Log.Message("Difference: %d (%f)", difference, diffPerc);

//		printf("%s\t%d\t%d\%d\t%d\t%d\t%d\t%d\t%f\t%d\t%f\n", group->fullRead->name,
//				group->fullRead->ReadId,
//				group->fullRead->length,
//				group->readNumber,
//				first, last,
//				distOnRead, distOnRef,
//				coveredOnRead,
//				difference, diffPerc);

//If difference < 0.1. assume that read doesn't span larger SVs and map
//using alignment
		if (pacbioDebug)
		Log.Message("%f < %f", fabs(diffPerc), 0.1f);
		if(fabs(diffPerc) < 0.1) {
			//Normal read. No event.
			MappedRead * read = group->fullRead;

			LocationScore * tmp = new LocationScore();

			if(startPos > endPos) {
				tmp->Location.m_Location = endPos;
				tmp->Location.setReverse(true);
			} else {
				tmp->Location.m_Location = startPos;
				tmp->Location.setReverse(false);
			}

			read->AllocScores(tmp, 1);
			read->Alignments = new Align[1];

			Timer tmr;
			tmr.ST();
			int corridor = std::max(abs(difference) * 2, (int)(read->length * 0.2));
			//corridor = corridor / 2 * 2;

			if(isReverse) {
				if (pacbioDebug)
				Log.Message("Read mapped reverse");

				read->computeReverseSeq();

				size_t const readSeqLen = endPosOnRead - startPosOnRead;
//				char * const readSeq = read->RevSeq + (read->length - endPosOnRead);
				char * const readSeq = new char[readSeqLen + 1];
				strncpy(readSeq, read->RevSeq + (read->length - endPosOnRead), readSeqLen);
				readSeq[readSeqLen] = '\0';
				if (pacbioDebug)
				Log.Message("ReadSeqLen: %d", readSeqLen);

				int const QEnd = startPosOnRead;
				int const QStart = read->length - endPosOnRead;

				//test-align python (ngila)
//				printf("%s\t%d\t", read->name, 1);
				read->Alignments[0] = computeAlignment(read->Scores[0].Location.m_Location, corridor, readSeq, readSeqLen, QStart, QEnd, read->length);
			} else {
//				char * const readSeq = read->Seq + startPosOnRead;
				size_t const readSeqLen = endPosOnRead - startPosOnRead;

				char * const readSeq = new char[readSeqLen + 1];
				strncpy(readSeq, read->Seq + startPosOnRead, readSeqLen);
				readSeq[readSeqLen] = '\0';
				if (pacbioDebug)
				Log.Message("ReadSeqLen: %d", readSeqLen);

				if (pacbioDebug)
				Log.Message("Start pos: %d, Length: %d", startPosOnRead, readSeqLen);

				int const QStart = startPosOnRead;
				int const QEnd = read->length - endPosOnRead;

				//test-align python (ngila)
//				printf("%s\t%d\t", read->name, 0);
				read->Alignments[0] = computeAlignment(read->Scores[0].Location.m_Location, corridor, readSeq, readSeqLen, QStart, QEnd, read->length);
			}
			if (pacbioDebug)
			Log.Message("CIGAR: %s", read->Alignments[0].pBuffer1);

			//test-align python (ngila)
//			printf("\t%d\n", read->Alignments[0].PositionOffset);
			read->Scores[0].Location.m_Location += read->Alignments[0].PositionOffset;//- (corridor >> 1); handled in computeAlingment

			if (pacbioDebug)
			Log.Message("Alignment took %fs", tmr.ET());

			read->Calculated = 1;

			WriteRead(read, true);

		} else {
			WriteRead(group->fullRead, false);
		}
	}
//	getchar();
}

void AlignmentBuffer::SaveRead(MappedRead * read, bool mapped) {
	//if (!argos) {
	WriteRead(read, mapped);
//	} else {
//		if (mapped) {
//			//Convert mapping position to RefId and position
//			for (int i = 0; i < read->Calculated; ++i) {
//				//TODO: fix for -n > 1
//				//Instead of setting mapped to false set score to 0 and don't print it in the end
//				mapped = SequenceProvider.convert(read->Scores[i].Location);
//			}
//		}
//		OutputReadBuffer::getInstance().addRead(read, mapped);
//		OutputReadBuffer::getInstance().getNextRead(m_Writer);
//	}
}

void AlignmentBuffer::WriteRead(MappedRead* read, bool mapped) {
	static int const topn = Config.GetInt("topn");
	if (mapped) {
		//Convert mapping position to RefId and position
		for (int i = 0; i < read->Calculated; ++i) {
			//TODO: fix for -n > 1
			//Instead of setting mapped to false set score to 0 and don't print it in the end
			mapped = SequenceProvider.convert(read->Scores[i].Location);
		}
	}
	if (read->Paired != 0) {
		if (topn == 1) {
			if (read->Paired->HasFlag(NGMNames::DeletionPending)) {
				if (read->hasCandidates() && read->Paired->hasCandidates()) {
					LocationScore * ls1 = &read->Scores[0];
					LocationScore * ls2 = &read->Paired->Scores[0];
					int distance =
							(ls2->Location.m_Location > ls1->Location.m_Location) ?
									ls2->Location.m_Location
											- ls1->Location.m_Location
											+ read->length :
									ls1->Location.m_Location
											- ls2->Location.m_Location
											+ read->Paired->length;

					//int distance = abs(read->TLS()->Location.m_Location - read->Paired->TLS()->Location.m_Location);

					pairInsertCount += 1;
					if (ls1->Location.getrefId() != ls2->Location.getrefId()
							|| distance < _NGM::sPairMinDistance
							|| distance > _NGM::sPairMaxDistance
							|| ls1->Location.isReverse()
									== ls2->Location.isReverse()) {
						//						Log.Message("%d != %d || %d < _%d || %d > %d || %d == %d", ls1->Location.getrefId() , ls2->Location.getrefId(), distance, _NGM::sPairMinDistance, distance, _NGM::sPairMaxDistance, ls1->Location.isReverse(), ls2->Location.isReverse());
						read->SetFlag(NGMNames::PairedFail);
						read->Paired->SetFlag(NGMNames::PairedFail);
						brokenPairs += 1;
					} else {
						pairInsertSum += distance;
					}
				}
				m_Writer->WritePair(read, 0, read->Paired, 0);
			}
		} else {
			Log.Error("TopN > 1 is currently not supported for paired end reads.");
			Fatal();
		}
	} else {
		m_Writer->WriteRead(read, mapped);
	}
	if (pairInsertCount % 1000 == 0) {
		NGM.Stats->validPairs = (pairInsertCount - brokenPairs) * 100.0f / pairInsertCount;
		//		NGM.Stats->insertSize = tSum * 1.0f / (tCount - brokenPairs);
	}
	NGM.GetReadProvider()->DisposeRead(read);
}
