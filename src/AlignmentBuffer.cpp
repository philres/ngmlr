#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>

#include "seqan/EndToEndAffine.h"
#include "OutputReadBuffer.h"
#include "Timing.h"

using std::vector;

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
		overallTime += tmr.ET();
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

	int maxTries = 1;

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

//		if (pacbioDebug) {
//			Log.Message("Ref: %s", refSeq);
//			Log.Message("Read: %s", readSeq);
//		}

		int mode = 0;
		cigarLength = aligner->SingleAlign(mode, corridor,
				(char const * const ) refSeq, (char const * const ) readSeq,
				align, clipping);

		int aligned = readLength - (align.QStart - clipping[0])
				- (align.QEnd - clipping[1]);

		if (aligned < 0.1f * readLength) {
			align.Score = -1.0f;
		}

		if (cigarLength == -1
				|| (aligned < 0.05f * readLength && maxTries > 0)) {
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
			if (cigarLength != -1) {
				maxTries -= 1;
				cigarLength = -1;
			}
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

void AlignmentBuffer::printDotPlotLine(int const id, char const * const name,
		int const onReadStart, int const onReadStop, loc const onRefStart,
		loc const onRefStop, float const score, bool const isReverse,
		int const type, int const status) {
	if (stdoutPrintDotPlot) {
		printf("%d\t%s\t%d\t%d\t%llu\t%llu\t%f\t%d\t%d\t%d\n", id, name,
				onReadStart, onReadStop, onRefStart, onRefStop, score,
				isReverse, type, status);
	}
}

void AlignmentBuffer::printDotPlotLine(int const id, char const * const name,
REAL const m, REAL const b, REAL const r, float const score,
		bool const isReverse, int const type, int const status) {
	if (stdoutPrintDotPlot) {
		printf("%d\t%s\t%.*f\t%.*f\t%.*f\t%f\t%f\t%d\t%d\t%d\n", id, name,
		DBL_DIG, m, DBL_DIG, b, DBL_DIG, r, 0.0f, score, isReverse, type,
				status);
	}
}

/* Finds longest increasing subsequence of reference positions in the
 * anchor list. Anchors are ordered by read position! */
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
				loc iRef =
						(anchors[i].isReverse) ?
								-1 * anchors[i].onRef : anchors[i].onRef;
				loc jRef =
						(anchors[j].isReverse) ?
								-1 * anchors[j].onRef : anchors[j].onRef;
				loc refDiff = llabs(iRef - jRef);

				loc readDiff = abs(anchors[i].onRead - anchors[j].onRead);

				loc diff = llabs(refDiff - readDiff);

//				Log.Message("cLIS: (rev %d, %d) %lld;%lld -  %lld - %lld = %lld; read %lld; diff %lld", anchors[j].isReverse, anchors[i].isReverse , anchors[i].onRef, anchors[j].onRef,  iRef, jRef, diff, readDiff, diff);

				if (DP[j] + 1 > DP[i] // New LIS must be longer then best so far
				&& anchors[j].isReverse == anchors[i].isReverse //Anchors have to be on the same strand to be joined!
						&& (diff < std::max(refDiff, readDiff) * 0.25f // The distance between the two anchors on the read and on the reference should not exceed a threshold. Adhoc threshold of 15% (INDEL sequencing error was to littel). 25% works better.
								|| (anchors[i].onRead == anchors[j].onRead
										&& llabs(refDiff) < readPartLength) // if anchors stem from the same read position accept ref positions that are in close proximity (to avoid additional segments caused be anchors that are only a few bp off)
						) && refDiff < readPartLength * 2.0f // Distance between anchors must be smaller than twice the anchor length. With this, gaps break the segment even if they are on the same "diagonal"
								) {
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

bool AlignmentBuffer::isSameDirection(Interval a, Interval b) {
	return a.isReverse == b.isReverse;
}

bool isIntervalInCorridor(REAL k, REAL d, REAL corridor,
		AlignmentBuffer::Interval testee) {

	bool inCorridor = false;
	REAL y = testee.onReadStart;
	loc upperRefStart = (loc) round((y - (d + corridor)) / k);
	loc lowerRefStart = (loc) round((y - (d - corridor)) / k);

	if (upperRefStart < lowerRefStart) {
		REAL tmp = upperRefStart;
		upperRefStart = lowerRefStart;
		lowerRefStart = tmp;
	}

	inCorridor = inCorridor
			|| (testee.onRefStart <= upperRefStart
					&& testee.onRefStart >= lowerRefStart);
//	Log.Message("Start: %llu <= %llu <= %llu", testee.onRefStart, lowerRefStart, upperRefStart);

	y = testee.onReadStop;
	loc upperRefStop = (loc) round((y - (d + corridor)) / k);
	loc lowerRefStop = (loc) round((y - (d - corridor)) / k);

	if (upperRefStop < lowerRefStop) {
		REAL tmp = upperRefStop;
		upperRefStop = lowerRefStop;
		lowerRefStop = tmp;
	}

	inCorridor = inCorridor
			&& (testee.onRefStop <= upperRefStop
					&& testee.onRefStop >= lowerRefStop);
//	Log.Message("Stop: %llu <= %llu <= %llu", testee.onRefStop, lowerRefStop, upperRefStop);

	//TODO: add max distance on read

	return inCorridor;
}
// Checks whether a is compatible (is located in the "corridor"
// of b.
bool AlignmentBuffer::isCompatible(Interval a, Interval b) {

	bool isCompatible = false;
	// Trade off (at the moment):
	// Bigger corridor: higher chance that alignment won't span event (score too low)
	// and a part of the reads is not aligned
	// Smaller corridor: reads are split, SVs are not visible in one alignment
	// TODO: after alignment add check if whole read was aligned. if not,
	// realign unanligned part
	REAL corridorSize = 2048;

	// Short intervals don't reliably define a corridor.
	if ((b.onReadStop - b.onReadStart) > readPartLength) {

		if (a.isReverse == b.isReverse) {
			// a and b are on the same strand
			isCompatible = isIntervalInCorridor(b.m, b.b, corridorSize, a);
		} else {
			if (pacbioDebug) {
				Log.Message("Intervals not on same strand");
			}
			// a and b are on different strands
			// Check if reads spans an inversion
			// If read spans an inversion, the
			// inverted part must be added to the segment
			// in order to keep consolidateSegements
			// to join the non inverted parts
			// through the inversion.

			// Switches the direction of the testee
			// Maybe not enough. Maybe switching the
			// direction of the tester is needed as well
			Interval tmp = a;
			tmp.onRefStart = a.onRefStop;
			tmp.onRefStop = a.onRefStart;
			isCompatible = isIntervalInCorridor(b.m, b.b, corridorSize, tmp);
		}

	}

	return isCompatible;

//	bool overlapsOnRead = false;
//	bool overlapsOnRef = false;
//
//	bool isInCorridor = false;
//
//	double k1 = (double) (a.onReadStop - a.onReadStart)
//			/ (a.onRefStop - a.onRefStart) * 1.0f;
//	double d1 = (double) (a.onRefStart * 1.0f * a.onReadStop
//			- a.onReadStart * 1.0f * a.onRefStop) / (a.onRefStart - a.onRefStop)
//			* 1.0f;
//
//	double k2 = (double) (b.onReadStop - b.onReadStart)
//			/ (double) (b.onRefStop - b.onRefStart) * 1.0f;
//	double d2 = (double) (b.onRefStart * 1.0f * b.onReadStop
//			- b.onReadStart * 1.0f * b.onRefStop)
//			/ (double) (b.onRefStart - b.onRefStop) * 1.0f;
//
////	Log.Message("a: %d, b: %d, diff: %d", d1, d2, d2 - d1);
//
//	loc a_y0 = -(d1 / k1);
//	loc b_y0 = -(d2 / k2);
//
//	if (pacbioDebug)
//		Log.Message("k1: %f, d1: %f", k1, d1);
//	if (pacbioDebug)
//		Log.Message("k2: %f, d2: %f", k2, d2);
//
//	if (pacbioDebug) {
//		Log.Message("Read start - a: %llu, b: %llu, diff: %llu", a_y0, b_y0, llabs(b_y0 - a_y0));
//	}
//
//	isInCorridor = llabs(b_y0 - a_y0) < 4096;
//
//	return !overlapsOnRead && !overlapsOnRef && isInCorridor;
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

AlignmentBuffer::Interval AlignmentBuffer::mergeIntervals(Interval a,
		Interval b) {
	if (a.onReadStart > b.onReadStart) {
		a.onReadStart = b.onReadStart;
		a.onRefStart = b.onRefStart;

	}
	if (a.onReadStop < b.onReadStop) {
		a.onReadStop = b.onReadStop;
		a.onRefStop = b.onRefStop;
	}
//	a.onReadStart = std::min(a.onReadStart, b.onReadStart);
//	a.onReadStop = std::max(a.onReadStop, b.onReadStop);
//	a.onRefStart = std::min(a.onRefStart, b.onRefStart);
//	a.onRefStop = std::max(a.onRefStop, b.onRefStop);
	a.score += b.score;
	return a;
}

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

/* Takes all anchors found for the read and transforms it in to intervals.
 * Intervals are continues candidate mapping for a read or parts of the
 * read. Contains repeatedly computing cLIS and reconciling of the found
 * intervals. */
AlignmentBuffer::Interval * AlignmentBuffer::findSubsequences(char * name,
		int id, Anchor * allFwdAnchors, int allFwdAnchorsLength,
		Anchor * allRevAnchors, int allRevAnchorsLength, int readParts,
		int readLenth, int & intervalsIndex) {

	if (pacbioDebug) {
		Log.Message("Finding LIS for read %d", id);
	}

	//TODO: remove fixed length
	int const maxMappedSegementCount = 100;
	MappedSegment * segments = new MappedSegment[maxMappedSegementCount];
	size_t segementsIndex = 0;

	//Sort by position on read. Probably not necessary!!
	std::sort(allFwdAnchors, allFwdAnchors + allFwdAnchorsLength,
			AlignmentBuffer::sortAnchorOnRead);

	int const maxcLISRunNumber = 8;
	int cLISRunNumber = 0;
	bool finished = false;
	while (!finished) {
		if (allFwdAnchorsLength > 0) {

			//Print remaining elements
			if (pacbioDebug) {
				std::cerr << "Elements: (i:onRead:onref:isReverse) ";
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

			if (pacbioDebug) {
				std::cerr << "cLIS" << cLISRunNumber << ":     ";
			}
			//Remove LIS from candidates
			int posInLIS = lisLength - 1;
			allRevAnchorsLength = 0;

			REAL * regX = new REAL[allFwdAnchorsLength];
			REAL * regY = new REAL[allFwdAnchorsLength];
			int pointNumber = 0;

			for (int i = 0; i < allFwdAnchorsLength; ++i) {
				if (posInLIS >= 0 && i == lis[posInLIS]) {

					int onRead = allFwdAnchors[lis[posInLIS]].onRead;

					//Print current LIS
					if (pacbioDebug) {
						std::cerr << lis[posInLIS] << ", ";
					}

					isReverse = allFwdAnchors[lis[posInLIS]].isReverse;
					intervalScore += allFwdAnchors[lis[posInLIS]].score;

					if(isReverse) {
						if(onRead < minOnRead) {
							minOnRead = onRead;
							minOnRef = allFwdAnchors[lis[posInLIS]].onRef + readPartLength;
						}

						if((onRead + readPartLength) > maxOnRead) {
							maxOnRead = onRead + readPartLength;
							maxOnRef = allFwdAnchors[lis[posInLIS]].onRef;
						}

						printDotPlotLine(id, name,
								onRead,
								onRead + readPartLength,
								allFwdAnchors[lis[posInLIS]].onRef + readPartLength,
								allFwdAnchors[lis[posInLIS]].onRef,
								allFwdAnchors[lis[posInLIS]].score,
								allFwdAnchors[lis[posInLIS]].isReverse,
								DP_TYPE_CLIS + cLISRunNumber, DP_STATUS_OK);

					} else {

						if(onRead < minOnRead) {
							minOnRead = onRead;
							minOnRef = allFwdAnchors[lis[posInLIS]].onRef;
						}

						if((onRead + readPartLength) > maxOnRead) {
							maxOnRead = onRead + readPartLength;
							maxOnRef = allFwdAnchors[lis[posInLIS]].onRef + readPartLength;
						}

						printDotPlotLine(id, name,
								onRead,
								onRead + readPartLength,
								allFwdAnchors[lis[posInLIS]].onRef,
								allFwdAnchors[lis[posInLIS]].onRef + readPartLength,
								allFwdAnchors[lis[posInLIS]].score,
								allFwdAnchors[lis[posInLIS]].isReverse,
								DP_TYPE_CLIS + cLISRunNumber, DP_STATUS_OK);

					}

					regY[pointNumber] = onRead;
					if(isReverse) {
						regX[pointNumber] = allFwdAnchors[lis[posInLIS]].onRef + readPartLength;
//						regX[pointNumber] = allFwdAnchors[lis[posInLIS]].onRef;
//						regY[pointNumber] = readLenth - onRead;
					} else {
						regX[pointNumber] = allFwdAnchors[lis[posInLIS]].onRef;
					}

					pointNumber += 1;

					//Remove from remaining elements
					posInLIS -= 1;
				} else {
					//Add to second array
					allRevAnchors[allRevAnchorsLength++] = allFwdAnchors[i];
				}
			}

			// Linear regression for segment
			REAL m,b,r;
			linreg(pointNumber,regX,regY,&m,&b,&r);
			delete[] regX; regX = 0;
			delete[] regY; regY = 0;
			if(pacbioDebug) {
				Log.Message("Regression: m=%.*f b=%.*f r=%.*f\n", DBL_DIG,m,DBL_DIG,b,DBL_DIG,r);
			}

			//Find intervall that covers all anchors best
			Interval interval;

			interval.isReverse = isReverse;
			interval.score = intervalScore;

			interval.onReadStart = minOnRead;
			interval.onReadStop = maxOnRead;

			interval.onRefStart = minOnRef;
			interval.onRefStop = maxOnRef;

			interval.m = m;
			interval.b = b;
			interval.r = r;

			if(interval.m != 0.0) {
				printDotPlotLine(id, name,
						interval.m, interval.b, r,
						interval.score,
						interval.isReverse,
						DP_TYPE_SEQMENTS_REG + cLISRunNumber, DP_STATUS_NOCOORDS);

//				printDotPlotLine(id, name,
//						-1 * interval.m, -1 * interval.b, r,
//						interval.score,
//						interval.isReverse,
//						DP_TYPE_SEQMENTS_REG + cLISRunNumber * 2 + 1, DP_STATUS_NOCOORDS);
			}

			printDotPlotLine(id, name,
					interval.onReadStart,
					interval.onReadStop,
					interval.onRefStart,
					interval.onRefStop,
					interval.score,
					interval.isReverse,
					DP_TYPE_SEQMENTS + cLISRunNumber, DP_STATUS_OK);

			//TODO: check chromosome borders

			if (pacbioDebug) {
				interval.print();
			}

			if(abs(interval.onReadStop - interval.onReadStart) > 0 &&
					llabs(interval.onRefStart - interval.onRefStop) > 0) {
				if(!constructMappedSegements(segments, interval, segementsIndex)) {
					Log.Message("Too many intervals found (>1000). Ignoring all others");
					lisLength = 0;
					delete[] lis;
//					return intervals;
					Fatal();
				}
			} else {
				if(pacbioDebug) {
					Log.Message("Interval too short.");
				}
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

			cLISRunNumber += 1;
			if (cLISRunNumber == maxcLISRunNumber) {
				// Max number of runs reached
				// TODO: find better stop criteria!
				finished = true;
			}
		} else {
			// If no unprocessed anchors are found anymore
			finished = true;
		}
	}

	intervalsIndex = 0;
	Interval * intervals = consolidateSegments(segments, segementsIndex, intervalsIndex);
	delete[] segments;
	segments = 0;

	// Extend all intervals to account for unmapped anchors at the end
	for(int i = 0; i < intervalsIndex; ++i) {
		Interval & interval = intervals[i];
		// If possible add 512 bp to start end end
		// min to avoid extending beyond the ends of the read
		int addStart = std::min(512, interval.onReadStart);
		int addEnd = std::min(512, readLenth - interval.onReadStop);

		// Try to match the slope of the interval when extending
		double lengthRatio = (interval.onReadStop - interval.onReadStart) * 1.0f / llabs(interval.onRefStop - interval.onRefStart) * 1.0f;

		if(interval.isReverse) {
			interval.onReadStart -= addStart;
			interval.onRefStart += (int) round(addStart / lengthRatio);

			interval.onReadStop += addEnd;
			interval.onRefStop -= (int) round(addEnd / lengthRatio);
		} else {
			interval.onReadStart -= addStart;
			interval.onRefStart -= (int) round(addStart / lengthRatio);

			interval.onReadStop += addEnd;
			interval.onRefStop += (int) round(addEnd / lengthRatio);

		}
	}

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

bool inline isInversion(float const nm) {
	float const minIdentity = 0.0;
	float const maxIdentity = 0.75;

	return nm > minIdentity && nm < maxIdentity;
}

void printLocation(char const * title, uloc position) {
	int len = 0;
	SequenceLocation loc;
	loc.m_Location = position;
	SequenceProvider.convert(loc);
	Log.Message("%s %s:%llu", title, SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location);
}

bool AlignmentBuffer::alignmentCheckForInversion(int const inversionLength,
		const int refCheckLength, SequenceLocation inversionCheckLocation,
		uloc inversionMidpointOnRead, const char* const readName,
		int inversionNumber, char* fullReadSeq) {

	bool inversionVerified = false;
	int const readCheckLength = 50;

	// TODO: don't allocate and delete every time. Can be done once
	// and reused
	const int readSeqLenght = readCheckLength * 2 + 10;
	char* readSeq = new char[readSeqLenght];
	char* revreadSeq = new char[readSeqLenght];
	memset(readSeq, '\0', readSeqLenght * sizeof(char));
	memset(revreadSeq, '\0', readSeqLenght * sizeof(char));
	const int refSeqLength = inversionLength + 2 * refCheckLength;
	char* refSeq = new char[refSeqLength + 10];
	memset(refSeq, '\0', (refSeqLength + 10) * sizeof(char));
	SequenceProvider.DecodeRefSequence(refSeq, 0, inversionCheckLocation.m_Location, refSeqLength);
	// Extract sequence from read
	strncpy(readSeq, fullReadSeq + inversionMidpointOnRead - readCheckLength,
			readCheckLength * 2);
	// Reverse complement
	computeReverseSeq(readSeq, revreadSeq, readCheckLength * 2);
	if (pacbioDebug) {
		Log.Message("Ref:     %s", refSeq);
		Log.Message("Read:    %s", readSeq);
		Log.Message("RevRead: %s", revreadSeq);
	}
	if (printInvCandidateFa) {
		printf(">%s_%d/1\n%s\n", readName, inversionNumber, refSeq);
		printf(">%s_%d/2\n%s\n", readName, inversionNumber, revreadSeq);
	}
	IAlignment* aligner = new EndToEndAffine();
	float scoreFwd = 0.0f;
	float scoreRev = 0.0f;
	static const float minScore = Config.GetFloat(MATCH_BONUS) * 1.0f
			* readCheckLength / 4.0f;
	aligner->SingleScore(10, 0, readSeq, refSeq, scoreFwd, 0);
	aligner->SingleScore(10, 0, revreadSeq, refSeq, scoreRev, 0);
	delete[] refSeq;
	refSeq = 0;
	delete[] readSeq;
	readSeq = 0;
	delete[] revreadSeq;
	revreadSeq = 0;
	delete aligner;
	aligner = 0;
	if (pacbioDebug) {
		SequenceProvider.convert(inversionCheckLocation);
		int len = 0;
		Log.Message("Inversion check at position %s:%llu - fwd: %f, rev: %f, min: %f", SequenceProvider.GetRefName(inversionCheckLocation.getrefId(), len), inversionCheckLocation.m_Location, scoreFwd, scoreRev, minScore);
	}
	inversionVerified = scoreRev > scoreFwd && scoreRev > minScore;
	//		inversionVerified = scoreRev > scoreFwd || scoreFwd < minScore;

	return inversionVerified;
}

bool AlignmentBuffer::validateInversion(Align const align,
		AlignmentBuffer::Interval const interval, int startInv, int stopInv,
		int startInvRead, int stopInvRead,
		char * fullReadSeq, AlignmentBuffer::Interval & leftOfInv,
		AlignmentBuffer::Interval & rightOfInv, MappedRead * read) {

	int len = 0;

	int inversionNumber = 0;

	int const minInversionLength = 25;
	// Length of regions that will be checked left and right to the midpoint of the
	// inversion on the reference
	int const refCheckLength = 250;

	int inversionLength = abs(stopInv - startInv);

	if (inversionLength > minInversionLength) {
		inversionNumber += 1;

		bool inversionVerified = false;

		uloc inversionMidpointOnRef = (startInv + stopInv) / 2;
		uloc inversionMidpointOnRead = (startInvRead + stopInvRead) / 2;

		SequenceLocation inversionCheckLocation;
		inversionCheckLocation.m_Location = interval.onRefStart
				+ align.PositionOffset + inversionMidpointOnRef
				- refCheckLength;

		inversionVerified = alignmentCheckForInversion(inversionLength,
				refCheckLength, inversionCheckLocation, inversionMidpointOnRead,
				read->name, inversionNumber, fullReadSeq);

		if (inversionVerified) {
			if (pacbioDebug) {
				Log.Message("Inversion detected: %llu, %llu", inversionMidpointOnRef,
						inversionMidpointOnRead);
			}

			if (interval.isReverse) {

				// alignStartOnRead is not relative to the full read, but to the part (interval)
				// that was aligned in the first place (the alignemnt that is scanned for an inversion)
				// Qstart is realtive to the full read.
				// However leftOfInv must be relative to the full read as well.
				// To define the stop position (start on the fwd strand) QStart can be used
				// For the start position (stop on the fwd strand) we use the position were
				// the inversion was detected. Since this position in relative to the aligned part
				// of the read (interval) we have to add additionalQstart
				int additionalQStart = align.QStart - align.firstPosition.readPosition;
				leftOfInv.onReadStop = read->length - (align.QStart);
				leftOfInv.onReadStart = read->length - (additionalQStart + inversionMidpointOnRead);

				// Position on the reference works perfectly with the mapping position
				// from the intial alignment. No need to account for the fact that only
				// a part of the read was aligned
				leftOfInv.onRefStart = interval.onRefStart + align.PositionOffset + align.firstPosition.refPosition;
				leftOfInv.onRefStop = interval.onRefStart + align.PositionOffset + inversionMidpointOnRef;

				leftOfInv.isReverse = interval.isReverse;

				// Same as for leftOfInv
				rightOfInv.onReadStart = read->length - (align.lastPosition.readPosition + additionalQStart);
				rightOfInv.onReadStop = read->length - (inversionMidpointOnRead + additionalQStart);

				rightOfInv.onRefStart = interval.onRefStart + align.PositionOffset + inversionMidpointOnRef;
				rightOfInv.onRefStop = interval.onRefStart + align.PositionOffset + align.lastPosition.refPosition;

				rightOfInv.isReverse = interval.isReverse;
			} else {
				leftOfInv.onReadStart = interval.onReadStart + align.firstPosition.readPosition;
				leftOfInv.onReadStop = interval.onReadStart + inversionMidpointOnRead;

				leftOfInv.onRefStart = interval.onRefStart + align.PositionOffset + align.firstPosition.refPosition;
				leftOfInv.onRefStop = interval.onRefStart + align.PositionOffset + inversionMidpointOnRef;
				leftOfInv.isReverse = interval.isReverse;

				rightOfInv.onReadStart = interval.onReadStart + inversionMidpointOnRead;
				rightOfInv.onReadStop = interval.onReadStart + align.lastPosition.readPosition;

				rightOfInv.onRefStart = interval.onRefStart + align.PositionOffset + inversionMidpointOnRef;
				rightOfInv.onRefStop = interval.onRefStart + align.PositionOffset + align.lastPosition.refPosition;
				rightOfInv.isReverse = interval.isReverse;
			}

			return true;
		}
	} else {
		if(pacbioDebug) {
			Log.Message("Inversion too short!");
		}
	}
	return false;
}

bool AlignmentBuffer::inversionDetection(Align const align,
		Interval const alignedInterval, char * readPartSeq,
		Interval & leftOfInv, Interval & rightOfInv, MappedRead * read) {

	static bool const enabled = Config.GetInt(NOINVERSIONS) != 1;
	if (!enabled) {
		return false;
	}

	//*********************//
	//Detect inversions
	//*********************//

	//Half of window size used to compute NM in SWCPU.cpp
	int const inversionPositionShift = 0;

	//Peaks (isInversion(nm) == true) that are less than maxDistance apart will be merged
	int const maxDistance = 20;
	int distance = maxDistance;

	int startInv = -1;
	int stopInv = -1;

	int startInvRead = -1;
	int stopInvRead = -1;

	SequenceLocation mappingLocation;
	mappingLocation.m_Location = alignedInterval.onRefStart + align.PositionOffset;
	SequenceProvider.convert(mappingLocation);

	if (stdoutErrorProfile) {
		int len = 0;
		for (int i = 0; i < align.alignmentLength; ++i) {
			printf("%s\t%llu\t%d\t%s\n", SequenceProvider.GetRefName(mappingLocation.getrefId(), len), mappingLocation.m_Location + align.nmPerPosition[i].refPosition, align.nmPerPosition[i].nm, read->name);
		}
	}

	// Extremely simple "peak-finding" algorithm
	int len = 0;
	for (int i = 0; i < align.alignmentLength; ++i) {
		float nm = (32 - align.nmPerPosition[i].nm) / 32.0f;

		if (startInv == -1) {
			if (isInversion(nm)) {
				startInv = align.nmPerPosition[i].refPosition
						- inversionPositionShift;
				startInvRead = align.nmPerPosition[i].readPosition
						- inversionPositionShift;
				stopInv = align.nmPerPosition[i].refPosition
						- inversionPositionShift;
				stopInvRead = align.nmPerPosition[i].readPosition
						- inversionPositionShift;
			}

		} else {	// if(stopInv == -1) {
			if (isInversion(nm)) {
				stopInv = align.nmPerPosition[i].refPosition
						- inversionPositionShift;
				stopInvRead = align.nmPerPosition[i].readPosition
						- inversionPositionShift;
				distance = maxDistance;
			} else {

				if (distance == 0) {

					if (stdoutInversionBed) {
						printf("%s\t%llu\t%llu\t%s\t%d\n",
						SequenceProvider.GetRefName(mappingLocation.getrefId(), len),
						mappingLocation.m_Location + startInv,
						mappingLocation.m_Location + stopInv, read->name, 0);
					}
					if (pacbioDebug) {
						Log.Message("Inversion detected: %d - %d, %d - %d (length: %d) on read %s",
						startInv, stopInv, startInvRead, stopInvRead,
						abs(stopInv - startInv), read->name);
					}
					if (validateInversion(align, alignedInterval, startInv, stopInv,
							startInvRead, stopInvRead, readPartSeq,
							leftOfInv, rightOfInv, read)) {
						return true;
					}

					startInv = -1;
					stopInv = -1;
					distance = maxDistance;
				} else {
					distance -= 1;
				}
			}
		}
	}

	return false;
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

	//Log.Message("%d %d", abs(interval.onReadStart - interval.onReadStop), llabs(interval.onRefStart - interval.onRefStop));
	if (llabs(interval.onReadStart - interval.onReadStop) == 0
			|| llabs(interval.onRefStart - interval.onRefStop) == 0) {
		Align align;
		return align;
	}

	Timer alignTimer;
	alignTimer.ST();
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

	alignTime += alignTimer.ET();
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
		Interval rightOfInv, MappedRead * read, Align * tmpAling,
		int & alignIndex, LocationScore * tmp, int mq) {

	if (pacbioDebug) {
		Log.Message("Inversion detected!");
		Log.Message("Interval:");
		interval.print();
	}
	/**********************************************/
	//Align part of read that is left of inversion
	/**********************************************/
	if (pacbioDebug) {
		Log.Message("Left Interval:");
		leftOfInv.print();
	}
	size_t readSeqLen = leftOfInv.onReadStop - leftOfInv.onReadStart;
	char* readSeq = extractReadSeq(readSeqLen, leftOfInv, read);
	Align align = alignInterval(read, leftOfInv, readSeq, readSeqLen);
	delete[] readSeq;
	readSeq = 0;
	delete[] align.nmPerPosition;
	align.nmPerPosition = 0;

	// Defines the inverted part of the alignment
	// start and end are computed from the alignments
	// left and right to it
	Interval inv;
	if (align.Score > 0.0f) {
		align.MQ = mq;
		tmpAling[alignIndex] = align;

		tmp[alignIndex].Location.m_Location = leftOfInv.onRefStart
		+ align.PositionOffset;	//- (corridor >> 1); handled in computeAlingment
		tmp[alignIndex].Location.setReverse(leftOfInv.isReverse);
		tmp[alignIndex].Score.f = align.Score;
		/**********************************************/
		//QEnd from first read part, gives length of
		//inversion that was covered by left part
		/**********************************************/
		// Set start of inverted read part
		inv.onReadStart = read->length - align.QEnd;
		inv.onRefStart = tmp[alignIndex].Location.m_Location
		+ align.lastPosition.refPosition;
		inv.isReverse = !leftOfInv.isReverse;

		read->Calculated += 1;
		alignIndex += 1;

	} else {
		if (pacbioDebug) {
			Log.Message("Alignment failed");
		}
		return false;
	}

	/**********************************************/
	//Align part of read that is right of inversion
	/**********************************************/
	if (pacbioDebug) {
		Log.Message("Right Interval:");
		rightOfInv.print();
	}
	readSeqLen = rightOfInv.onReadStop - rightOfInv.onReadStart;
	readSeq = extractReadSeq(readSeqLen, rightOfInv, read);
	align = alignInterval(read, rightOfInv, readSeq, readSeqLen);
	delete[] readSeq;
	readSeq = 0;
	delete[] align.nmPerPosition;
	align.nmPerPosition = 0;

	if (align.Score > 0.0f) {
		align.MQ = mq;
		tmpAling[alignIndex] = align;

		tmp[alignIndex].Location.m_Location = rightOfInv.onRefStart
		+ align.PositionOffset;	//- (corridor >> 1); handled in computeAlingment
		tmp[alignIndex].Location.setReverse(rightOfInv.isReverse);
		tmp[alignIndex].Score.f = align.Score;
		// Set end of inverted read part
		/**********************************************/
		//QStart from second read part, gives length
		//of inversion that was covered by second read
		/**********************************************/
		inv.onReadStop = align.QStart;
		inv.onRefStop = tmp[alignIndex].Location.m_Location
		+ align.firstPosition.refPosition;

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
	if(!inv.isReverse) {
		// Reverse coordinates on read
		uloc tmp = read->length - inv.onReadStart;
		inv.onReadStart = read->length - inv.onReadStop;
		inv.onReadStop = tmp;
	}
	// Extend
	int extLength = (int) round((inv.onReadStop - inv.onReadStart) * 0.2f);
	inv.onReadStart -= extLength;
	inv.onReadStop += extLength;
	inv.onRefStart -= extLength;
	inv.onRefStop += extLength;

	int inversionLength = inv.onReadStop - inv.onReadStart;
	if (pacbioDebug) {
		Log.Message("Lenght of inversion is %d", inversionLength);
		Log.Message("Inverted Interval:");
		inv.print();
	}
	if(inversionLength > 40) {
		readSeqLen = inv.onReadStop - inv.onReadStart;
		readSeq = extractReadSeq(readSeqLen, inv, read);
		align = alignInterval(read, inv, readSeq, readSeqLen);
		delete[] readSeq;
		readSeq = 0;
		delete[] align.nmPerPosition;
		align.nmPerPosition = 0;

		if (align.Score > 0.0f) {
			align.MQ = mq;
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
	} else {
		Log.Message("Inversion skipped");
	}

	return true;
}

void AlignmentBuffer::alignSingleOrMultipleIntervals(MappedRead * read,
		Interval interval, LocationScore * tmp, Align * tmpAling,
		int & alignIndex) {

	if (interval.onRefStart > interval.onRefStop) {
		loc tmp = interval.onRefStart;
		interval.onRefStart = interval.onRefStop;
		interval.onRefStop = tmp;
	}

	size_t const readSeqLen = interval.onReadStop - interval.onReadStart;
	char * readPartSeq = extractReadSeq(readSeqLen, interval, read);

	Align align = alignInterval(read, interval, readPartSeq, readSeqLen);

	if (pacbioDebug) {
		SequenceLocation seqLoc2;
		seqLoc2.m_Location = interval.onRefStart + align.PositionOffset;
		SequenceProvider.convert(seqLoc2);
		Log.Message("Mapping Pos (%s): %llu -- %llu", read->name, interval.onRefStart + align.PositionOffset, seqLoc2.m_Location);
	}

	Interval leftOfInv;
	Interval rightOfInv;
	bool inversionDetected = false;
	bool inversionAligned = false;
	if ((inversionDetected = inversionDetection(align, interval, readPartSeq,
			leftOfInv, rightOfInv, read))) {

		int mq = computeMappingQuality(align, read->length);
		inversionAligned = alignInversion(interval, leftOfInv, rightOfInv, read,
				tmpAling, alignIndex, tmp, mq);
		if (!inversionAligned) {
			SequenceLocation loc;
			loc.m_Location = leftOfInv.onRefStop;
			SequenceProvider.convert(loc);
			SequenceLocation loc2;
			loc2.m_Location = rightOfInv.onRefStart;
			SequenceProvider.convert(loc2);
			int len = 0;
			Log.Message("Potential inversion detected at %s:%llu-%llu but not taken into account due to errors. (read: %s)", SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location, loc2.m_Location, read->name);
		}
	} else {
		if (pacbioDebug) {
			Log.Message("No inversion detected!");
		}
	}
	delete[] align.nmPerPosition;
	align.nmPerPosition = 0;

	if (!inversionDetected || !inversionAligned) {
		/**********************************************/
		// No inversion detected
		/**********************************************/
		if (align.Score > 0.0f) {
			align.MQ = computeMappingQuality(align, read->length);
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
//	TODO: use nmPosition instead (see IAlignment.h)
//		int QStart = ((int*) (align.ExtendedData))[1];
//		int QEnd = readSeqLen - ((int*) (align.ExtendedData))[3];
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

	delete[] readPartSeq;
	readPartSeq = 0;

}

int AlignmentBuffer::computeMappingQuality(Align const & alignment,
		int readLength) {
	std::vector<IntervalTree::Interval<int> > results;

	readCoordsTree->findOverlapping(alignment.QStart,
			readLength - alignment.QEnd, results);
//	Log.Message("On read: %d - %d", alignment.QStart, readLength - alignment.QEnd);
	int mqSum = 0;
	int mqCount = 0;
	for (int j = 0; j < results.size(); ++j) {
//		Log.Message("\tOverlap: %d;%d;%d", results[j].start, results[j].stop, results[j].value);
		mqSum += results[j].value;
		mqCount += 1;
	}
	return (int) (mqSum * 1.0f / mqCount);
}

//bool sortAlignPerRefPos(Align a,
//		Align b) {
//	return a. < b.onReadStart;
//}

bool sortMappedSegements(IntervalTree::Interval<AlignmentBuffer::Interval *> a,
		IntervalTree::Interval<AlignmentBuffer::Interval *> b) {
//	return a.value.onReadStart < b.value.onReadStart;
//	return a.value->score > b.value->score;
	return (a.value->onReadStop - a.value->onReadStart)
			> (b.value->onReadStop - b.value->onReadStart);
}

bool isFullyContainedOnRead(AlignmentBuffer::Interval * shortInterval,
		AlignmentBuffer::Interval * longInterval) {

//	int threshold = 0.01f * (longInterval->onReadStop - longInterval->onReadStart);
	int threshold = 10;
	return (longInterval->onReadStart - shortInterval->onReadStart) <= threshold
			&& (longInterval->onReadStop - shortInterval->onReadStop)
					>= -threshold;
}

bool isFullyContainedOnRef(AlignmentBuffer::Interval * shortInterval,
		AlignmentBuffer::Interval * longInterval) {

//	int threshold = 0.10f * (longInterval->onReadStop - longInterval->onReadStart);
	int threshold = 10;
	return (longInterval->onRefStart - shortInterval->onRefStart) <= threshold
			&& (longInterval->onRefStop - shortInterval->onRefStop)
					>= -threshold;
}

bool isFullyContained(AlignmentBuffer::Interval * shortInterval,
		AlignmentBuffer::Interval * longInterval) {

	return isFullyContainedOnRead(shortInterval, longInterval)
			&& isFullyContainedOnRef(shortInterval, longInterval);
}

bool isValidOverlap(AlignmentBuffer::Interval * shortInterval,
		AlignmentBuffer::Interval * longInterval) {
	return true;
}

bool isValidOverlapRef(AlignmentBuffer::Interval * shortInterval,
		AlignmentBuffer::Interval * longInterval) {
	return shortInterval->isReverse == longInterval->isReverse
			&& shortInterval->onRefStart <= longInterval->onRefStop
			&& shortInterval->onRefStop >= longInterval->onRefStart;
}

void AlignmentBuffer::reconcileRead(ReadGroup * group) {

	std::vector<IntervalTree::Interval<Interval *> > mappedSegements;

	MappedRead * read = group->fullRead;

	// Not the best strategy, however only used for debug
	// plots at the moment (ngm itself doesn't rely on read strand information)
	bool readIsReverse = read->Scores[0].Location.isReverse();

	if (pacbioDebug) {
		Log.Message("Reconciling read: %s (%d)", read->name, read->ReadId);
	}
	// Convert Score + Alignment info to Interval (mappedSegement) object
	// mappedSegemnts will be used to detect overlaps in alignments
	// TODO: create mapped segments in an early stage of the algorithm
	// maybe merge Align and Interval classes
	for (int i = 0; i < read->Calculated; ++i) {
		int diffOnRef = read->Alignments[i].lastPosition.refPosition
				- read->Alignments[i].firstPosition.refPosition;

		Interval * mappedSegement = new Interval();
		mappedSegement->id = i;
		mappedSegement->onRefStart = read->Scores[i].Location.m_Location;
		mappedSegement->onRefStop = read->Scores[i].Location.m_Location
				+ diffOnRef;
		mappedSegement->isReverse = read->Scores[i].Location.isReverse();
		mappedSegement->isProcessed = false;
		mappedSegement->score = read->Alignments[i].Score;
		if (mappedSegement->isReverse) {
			// If mapped to the reverse strand, QStart and QEnd have to be translated to
			// plus strand
			mappedSegement->onReadStart = read->Alignments[i].QEnd;
			mappedSegement->onReadStop = read->length
					- read->Alignments[i].QStart - 1;
		} else {
			mappedSegement->onReadStart = read->Alignments[i].QStart;
			mappedSegement->onReadStop = read->length - read->Alignments[i].QEnd
					- 1;
		}
		mappedSegements.push_back(
				IntervalTree::Interval<Interval *>(mappedSegement->onReadStart,
						mappedSegement->onReadStop, mappedSegement));

		if (stdoutPrintMappedSegments) {
			fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", read->name, i,
					read->length, mappedSegement->onReadStart,
					mappedSegement->onReadStop,
					read->Scores[i].Location.isReverse(),
					read->Alignments[i].MQ, read->Alignments[i].Score);
		}
	}

	std::sort(mappedSegements.begin(), mappedSegements.end(),
			sortMappedSegements);
	for (int i = 0; i < mappedSegements.size(); ++i) {
		if (pacbioDebug) {
			Log.Message("%d: Read: %d - %d, Ref: %llu - %llu (%d) - %f - %d", mappedSegements[i].value->id, mappedSegements[i].value->onReadStart, mappedSegements[i].value->onReadStop, mappedSegements[i].value->onRefStart, mappedSegements[i].value->onRefStop, mappedSegements[i].value->isReverse, mappedSegements[i].value->score, mappedSegements[i].value->isProcessed);
		}

		if(!readIsReverse) {
			printDotPlotLine(read->ReadId, read->name,
					mappedSegements[i].value->onReadStart,
					mappedSegements[i].value->onReadStop,
					mappedSegements[i].value->onRefStart,
					mappedSegements[i].value->onRefStop,
					mappedSegements[i].value->score,
					mappedSegements[i].value->isReverse,
					DP_TYPE_RESULT + mappedSegements[i].value->id, DP_STATUS_OK);
		} else {
			printDotPlotLine(read->ReadId, read->name,
					mappedSegements[i].value->onReadStop,
					mappedSegements[i].value->onReadStart,
					mappedSegements[i].value->onRefStart,
					mappedSegements[i].value->onRefStop,
					mappedSegements[i].value->score,
					mappedSegements[i].value->isReverse,
					DP_TYPE_RESULT + mappedSegements[i].value->id, DP_STATUS_OK);
		}
	}

	IntervalTree::IntervalTree<Interval *> mappedSegmentsTree(mappedSegements);

	bool overlapFound = false;
	for (int i = 0; i < mappedSegements.size(); ++i) {
		if (!mappedSegements[i].value->isProcessed) {
			std::vector<IntervalTree::Interval<Interval *> > results;
			mappedSegements[i].value->isProcessed = true;
			mappedSegmentsTree.findOverlapping(
					mappedSegements[i].value->onReadStart,
					mappedSegements[i].value->onReadStop, results);
			for (int j = 0; j < results.size(); ++j) {
				if (/*!(results[j].value == mappedSegements[i].value) &&*/!results[j].value->isProcessed
						&& isValidOverlap(results[j].value,
								mappedSegements[i].value)) {
					if (isFullyContainedOnRead(results[j].value,
							mappedSegements[i].value)) {
						if (isFullyContainedOnRef(results[j].value,
								mappedSegements[i].value)) {
							if (pacbioDebug) {
								Log.Message("%d is fully contained in %d on read and ref", results[j].value->id, mappedSegements[i].value->id);
							}
							read->Alignments[results[j].value->id].skip = true;
							results[j].value->isProcessed = true;
						} else {
							read->Alignments[results[j].value->id].MQ = 0;
							if(pacbioDebug) {
								Log.Message("%d is fully contained in %d on read", results[j].value->id, mappedSegements[i].value->id);
							}
							results[j].value->isProcessed = true;
						}
					} else {
						if(isValidOverlapRef(results[j].value,
										mappedSegements[i].value)) {
							Log.Error("Mapped segment %d overlaps with %d on read %s and reference", mappedSegements[i].value->id, results[j].value->id, read->name);
							//Fatal();
							results[j].value->isProcessed = true;
							overlapFound = true;

						}
					}
				}
			}
		}
	}

	if (overlapFound) {
//		getchar();
	}

	for (int i = 0; i < read->Calculated; ++i) {
		if (!read->Alignments[mappedSegements[i].value->id].skip) {
			if (!readIsReverse) {
				printDotPlotLine(read->ReadId, read->name,
						mappedSegements[i].value->onReadStart,
						mappedSegements[i].value->onReadStop,
						mappedSegements[i].value->onRefStart,
						mappedSegements[i].value->onRefStop,
						mappedSegements[i].value->score,
						mappedSegements[i].value->isReverse,
						DP_TYPE_RESULT_CONS + mappedSegements[i].value->id,
						DP_STATUS_OK);
			} else {
				printDotPlotLine(read->ReadId, read->name,
						mappedSegements[i].value->onReadStop,
						mappedSegements[i].value->onReadStart,
						mappedSegements[i].value->onRefStart,
						mappedSegements[i].value->onRefStop,
						mappedSegements[i].value->score,
						mappedSegements[i].value->isReverse,
						DP_TYPE_RESULT_CONS + mappedSegements[i].value->id,
						DP_STATUS_OK);
			}
		}
	}
}

void AlignmentBuffer::processLongReadLIS(ReadGroup * group) {

	Timer tmr;
	tmr.ST();

	// Parts of read that were aligned plus MQ of subalignments
	std::vector<IntervalTree::Interval<int> > treeIntervals;

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

	float fwdScore = 0.0f;
	float revScore = 0.0f;
	for (int j = 0; j < group->readNumber; ++j) {
		MappedRead * part = group->reads[j];

		int onRead = j;

		treeIntervals.push_back(
				IntervalTree::Interval<int>(onRead * readPartLength,
						onRead * readPartLength + readPartLength,
						part->mappingQlty));

		float minScore =
				(part->numScores() > 0) ? part->Scores[0].Score.f * 0.8 : 0.0f;
		minScore = std::max(minScore, minGroupScore);
//		minScore = 0.0f;

		// Get all mapping positions from anchors (non overlapping 512bp parts of reads)
		// and convert them to Anchor objects
		// TODO: remove fixed length of 100
		int maxCandidatesPerReadPart = 100;
		for (int k = 0; k < part->numScores(); ++k) {

			// If anchor has to many mapping positions -> ignore
			if (part->numScores() < maxCandidatesPerReadPart) {
				if (part->Scores[k].Score.f > minScore) {
					// Anchor is valid and will be used
					Anchor & anchor = anchorsFwd[anchorFwdIndex++];

					anchor.score = part->Scores[k].Score.f;
					anchor.isReverse = part->Scores[k].Location.isReverse();
					if (anchor.isReverse) {
						revScore += anchor.score;
					} else {
						fwdScore += anchor.score;
					}
					anchor.type = DP_STATUS_OK;

					// It would be best to convert reads or read parts that map to the negative strand
					// Immediately to plus strand.
					// If somebody tells me how to do this properly, I'll happily change this.
					// For now Anchors can be an the plus or minus strand!
					// Problem: Transforming coordinates from negative anchors to plus strand
					// is not possible without knowing which part of the read is from plus
					// and which form minus strand. At this point it would be hard to infer that
					// (although possible)
					if (anchor.isReverse) {
						anchor.onRead = onRead * readPartLength;
						//anchor.onRead = group->fullRead->length - (onRead + 1) * readPartLength;
						anchor.onRef = part->Scores[k].Location.m_Location;
						printDotPlotLine(group->fullRead->ReadId,
								group->fullRead->name, anchor.onRead,
								anchor.onRead + readPartLength,
								part->Scores[k].Location.m_Location
										+ readPartLength,
								part->Scores[k].Location.m_Location,
								part->Scores[k].Score.f,
								part->Scores[k].Location.isReverse(),
								DP_TYPE_UNFILTERED, DP_STATUS_OK);
					} else {
						anchor.onRead = onRead * readPartLength;
						//anchor.onRead = group->fullRead->length - (onRead + 1) * readPartLength;
						anchor.onRef = part->Scores[k].Location.m_Location;
						printDotPlotLine(group->fullRead->ReadId,
								group->fullRead->name, anchor.onRead,
								anchor.onRead + readPartLength,
								part->Scores[k].Location.m_Location,
								part->Scores[k].Location.m_Location
										+ readPartLength,
								part->Scores[k].Score.f,
								part->Scores[k].Location.isReverse(),
								DP_TYPE_UNFILTERED, DP_STATUS_OK);
					}

					if (pacbioDebug) {
						Log.Message("\t%d\t%f at %llu", j, part->Scores[k].Score.f, part->Scores[k].Location.m_Location);
					}

				} else {
					//Score too low
					// Still print for visualization
					printDotPlotLine(group->fullRead->ReadId,
							group->fullRead->name, onRead, onRead + readPartLength,
							part->Scores[k].Location.m_Location,
							part->Scores[k].Location.m_Location + readPartLength,
							part->Scores[k].Score.f,
							part->Scores[k].Location.isReverse(),
							DP_TYPE_UNFILTERED, DP_STATUS_LOWSCORE);
				}
			} else {
				// Repetitive
				// Still print for visualization
				printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
						onRead, onRead + readPartLength,
						part->Scores[k].Location.m_Location,
						part->Scores[k].Location.m_Location + readPartLength,
						part->Scores[k].Score.f,
						part->Scores[k].Location.isReverse(),
						DP_TYPE_UNFILTERED,
						DP_STATUS_REPETITIVE);
			}
		}

		if (part->numScores() == 0) {
			//No hits found
			if (pacbioDebug) {
				Log.Message("No hits found for part %d", onRead);
			}
			printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
					onRead, onRead + readPartLength, 0, 0, 0.0f, 0, DP_TYPE_UNFILTERED, DP_STATUS_NOHIT);
		}
	}

	// This is not used at the moment. Reads are mapped indipendent of their strand
	if (pacbioDebug) {
		Log.Message("Trying to determine read strand: %f fwd, %f rev", fwdScore, revScore);
		if (revScore > fwdScore) {
			Log.Message("Read is 'probably' from reverse strand");
		} else {
			Log.Message("Read is 'probably' from forward strand");
		}
	}

	readCoordsTree = new IntervalTree::IntervalTree<int>(treeIntervals);

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
			printDotPlotLine(read->ReadId, read->name,
					intervals[i].onReadStart,
					intervals[i].onReadStop,
					intervals[i].onRefStart,
					intervals[i].onRefStop,
					intervals[i].score,
					intervals[i].isReverse,
					DP_TYPE_SEQMENTS_CONS + i, DP_STATUS_OK);

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

	delete readCoordsTree;
	readCoordsTree = 0;

	processTime += tmr.ET();
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

//////////////////////////////////////////////////
///////////////// Deprecated /////////////////////
//////////////////////////////////////////////////

//AlignmentBuffer::Interval AlignmentBuffer::splitInterval(Interval a,
//		Interval b) {
//
//	Log.Message("=======================");
//	a.print();
//	b.print();
////	getchar();
//	return a;
//}

//void AlignmentBuffer::processLongRead(ReadGroup * group) {
////					Log.Message("Read group with id %d finished", group->readId);
////					Log.Message("Name: %s", cur_read->name);
////					Log.Message("Reads in group: %d", group->readNumber);
////					Log.Message("Reads finished: %d", group->readsFinished);
////					Log.Message("Fwd: %d, Rev: %d", group->fwdMapped, group->reverseMapped);
////					Log.Message("Avg best score %f", group->bestScoreSum * 1.0f / group->readsFinished);
//
//	Log.Error("Should not be here. AlignmentBuffer processLongRead");
//	Fatal();
//
//	float avgGroupScore = group->bestScoreSum * 1.0f / group->readsFinished;
//	float minGroupScore = avgGroupScore * 0.8f;
//
////	for (int j = 0; j < group->readNumber; ++j) {
////		MappedRead * part = group->reads[j];
////		//Log.Message("ID: %d (has %d scores)", part->ReadId, part->numScores());
////		float minScore = part->Scores[0].Score.f * 0.8;
////		for (int k = 0; k < part->numScores(); ++k) {
////			if (part->Scores[k].Score.f > minScore) {
////				//Log.Message("\t%f at %llu", part->Scores[k].Score.f, part->Scores[k].Location.m_Location);
////			}
////		}
////	}
//
////Find first read part that maps with a min score
//	int first = 0;
//	while (first < group->readNumber
//			&& (group->reads[first]->numScores() == 0
//					|| group->reads[first]->Scores[0].Score.f < minGroupScore)) {
//		first += 1;
//	}
//
////Find last read part that maps with a min score
//	int last = group->readNumber - 1;
//	while (last >= 0
//			&& (group->reads[last]->numScores() == 0
//					|| group->reads[last]->Scores[0].Score.f < minGroupScore)) {
//		last -= 1;
//	}
//
//	if (first == group->readNumber || last < 0 || group->readNumber < 3) {
//		Log.Verbose("Could not map read.");
//	} else {
////Distance on read between start of first and last mapped read part
////+1 to take the full last read part into account
//		int distOnRead = (last - first + 1) * 512;
//
////If not the whole read is aligned add half of the read part size to alignment
//		size_t startPosOnRead = std::max(0, first * 512 - 256);
//		size_t endPosOnRead = std::min(group->fullRead->length, (last + 1) * 512 + 256);
//
//		if (pacbioDebug)
//		Log.Message("Name: %s", group->fullRead->name);
//		if (pacbioDebug)
//		Log.Message("On read (length %d): %d to %d (dist %d)", group->fullRead->length, startPosOnRead, endPosOnRead, distOnRead);
//
//		bool isReverse = false;
////Compute distance between mapped location of first and last mapped read part on reference
//		uloc endPos = group->reads[last]->Scores[0].Location.m_Location;
//		uloc startPos = group->reads[first]->Scores[0].Location.m_Location;
//		int distOnRef = 0;
//		if(startPos > endPos) {
//			distOnRef = startPos - endPos;
//			isReverse = true;
//		} else {
//			distOnRef = endPos - startPos;
//		}
////+512 to take full length of last read part into account
//		distOnRef += 512;
//		if (pacbioDebug)
//		Log.Message("Start pos on ref: %llu to %llu (dist %llu)", startPos, endPos, distOnRef);
//
////Log.Message("Start: %llu, End: %llu", first, last);
////Log.Message("Read length: %d", group->fullRead->length);
//
//		float coveredOnRead = (distOnRead) * 100.0f / group->fullRead->length;
//		if (pacbioDebug)
//		Log.Message("Covered on read: %f", coveredOnRead);
////Log.Message("On read: %d, On ref: %llu", distOnRead, distOnRef);
//
////Difference between distance in read cooridnates and distance in ref coordinates
////If read doesn't span larger structural variations, difference should only be
////caused by PacBio sequence error model
//		int difference = distOnRead - distOnRef;
//		float diffPerc = (distOnRead - distOnRef) * 1.0f / group->fullRead->length;
//		if (pacbioDebug)
//		Log.Message("Difference: %d (%f)", difference, diffPerc);
//
////		printf("%s\t%d\t%d\%d\t%d\t%d\t%d\t%d\t%f\t%d\t%f\n", group->fullRead->name,
////				group->fullRead->ReadId,
////				group->fullRead->length,
////				group->readNumber,
////				first, last,
////				distOnRead, distOnRef,
////				coveredOnRead,
////				difference, diffPerc);
//
////If difference < 0.1. assume that read doesn't span larger SVs and map
////using alignment
//		if (pacbioDebug)
//		Log.Message("%f < %f", fabs(diffPerc), 0.1f);
//		if(fabs(diffPerc) < 0.1) {
//			//Normal read. No event.
//			MappedRead * read = group->fullRead;
//
//			LocationScore * tmp = new LocationScore();
//
//			if(startPos > endPos) {
//				tmp->Location.m_Location = endPos;
//				tmp->Location.setReverse(true);
//			} else {
//				tmp->Location.m_Location = startPos;
//				tmp->Location.setReverse(false);
//			}
//
//			read->AllocScores(tmp, 1);
//			read->Alignments = new Align[1];
//
//			Timer tmr;
//			tmr.ST();
//			int corridor = std::max(abs(difference) * 2, (int)(read->length * 0.2));
//			//corridor = corridor / 2 * 2;
//
//			if(isReverse) {
//				if (pacbioDebug)
//				Log.Message("Read mapped reverse");
//
//				read->computeReverseSeq();
//
//				size_t const readSeqLen = endPosOnRead - startPosOnRead;
////				char * const readSeq = read->RevSeq + (read->length - endPosOnRead);
//				char * const readSeq = new char[readSeqLen + 1];
//				strncpy(readSeq, read->RevSeq + (read->length - endPosOnRead), readSeqLen);
//				readSeq[readSeqLen] = '\0';
//				if (pacbioDebug)
//				Log.Message("ReadSeqLen: %d", readSeqLen);
//
//				int const QEnd = startPosOnRead;
//				int const QStart = read->length - endPosOnRead;
//
//				//test-align python (ngila)
////				printf("%s\t%d\t", read->name, 1);
//				read->Alignments[0] = computeAlignment(read->Scores[0].Location.m_Location, corridor, readSeq, readSeqLen, QStart, QEnd, read->length);
//			} else {
////				char * const readSeq = read->Seq + startPosOnRead;
//				size_t const readSeqLen = endPosOnRead - startPosOnRead;
//
//				char * const readSeq = new char[readSeqLen + 1];
//				strncpy(readSeq, read->Seq + startPosOnRead, readSeqLen);
//				readSeq[readSeqLen] = '\0';
//				if (pacbioDebug)
//				Log.Message("ReadSeqLen: %d", readSeqLen);
//
//				if (pacbioDebug)
//				Log.Message("Start pos: %d, Length: %d", startPosOnRead, readSeqLen);
//
//				int const QStart = startPosOnRead;
//				int const QEnd = read->length - endPosOnRead;
//
//				//test-align python (ngila)
////				printf("%s\t%d\t", read->name, 0);
//				read->Alignments[0] = computeAlignment(read->Scores[0].Location.m_Location, corridor, readSeq, readSeqLen, QStart, QEnd, read->length);
//			}
//			if (pacbioDebug)
//			Log.Message("CIGAR: %s", read->Alignments[0].pBuffer1);
//
//			//test-align python (ngila)
////			printf("\t%d\n", read->Alignments[0].PositionOffset);
//			read->Scores[0].Location.m_Location += read->Alignments[0].PositionOffset;//- (corridor >> 1); handled in computeAlingment
//
//			if (pacbioDebug)
//			Log.Message("Alignment took %fs", tmr.ET());
//
//			read->Calculated = 1;
//
//			WriteRead(read, true);
//
//		} else {
//			WriteRead(group->fullRead, false);
//		}
//	}
////	getchar();
//}

//void AlignmentBuffer::getSeqAnLIS(int allFwdAnchorsLength, int id,
//		Anchor* allFwdAnchors, char* name) {
//
//	//TODO: try longest common subsequence with 1 specific symbol per read part.
//	//Second sequence consists of all matching positions on the reference (orderd)
//	//Positions get the same symbol as a read part if they got a alignment score
//	String<loc> seqFwd;
//	String<loc> seqRev;
//	String<float> weightsFwd;
//	String<float> weightsRev;
//	String<loc, Block<> > posFwd;
//	String<loc, Block<> > posRev;
//	for (int i = 0; i < allFwdAnchorsLength; ++i) {
//		if (allFwdAnchors[i].isReverse) {
//			appendValue(seqFwd, allFwdAnchors[i].onRef * -1);
//		} else {
//			appendValue(seqFwd, allFwdAnchors[i].onRef);
//		}
//		appendValue(weightsFwd, allFwdAnchors[i].score);
//	}
//	//	for (int i = 0; i < allRevAnchorsLength; ++i) {
//	//		appendValue(seqRev, allRevAnchors[i].onRef);
//	//		appendValue(weightsRev, allRevAnchors[i].score);
//	//
//	//	}
//	longestIncreasingSubsequence(seqFwd, posFwd);
////	float w = heaviestIncreasingSubsequence(seqFwd, weightsFwd, posFwd);
//
//	if (pacbioDebug)
//		std::cerr << "SeqAN LIS: ";
//	for (int i = length(posFwd) - 1; i >= 0; --i) {
//		if (pacbioDebug)
//			std::cerr << posFwd[i] << ", ";
////		printDotPlotLine(id, name, allFwdAnchors[posFwd[i]].onRead * 512,
////				allFwdAnchors[posFwd[i]].onRead * 512 + 512,
////				allFwdAnchors[posFwd[i]].onRef,
////				allFwdAnchors[posFwd[i]].onRef + 512,
////				allFwdAnchors[posFwd[i]].score,
////				allFwdAnchors[posFwd[i]].isReverse, TYPE_LIS_SEQAN, STATUS_OK);
//	}
//	if (pacbioDebug)
//		std::cerr << std::endl;
//
//}

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

//bool AlignmentBuffer::inversionDetectionArndt(Align const align,
//		Interval const interval, int const readLength, char * fullReadSeq,
//		Interval & leftOfInv, Interval & rightOfInv, Interval & inv,
//		char const * const readName) {
//
//	static bool const enabled = Config.GetInt(NOINVERSIONS) != 1;
//
//	if (!enabled) {
//		return false;
//	}
//
//	loc start = interval.onRefStart;
//
//	SequenceLocation seqLoc;
//	seqLoc.m_Location = start + align.PositionOffset;
//	SequenceProvider.convert(seqLoc);
//
//	if (pacbioDebug) {
//
//		Log.Message("Mapping Pos: %llu", seqLoc.m_Location);
//	}
//
//	bool inversionFound = false;
//
////*********************//
////Detect inversions
////*********************//
//
//	//Not bp but differences (mismatch, insertion, deletion)
//	//Inversion is only detected if NM is above threshold for 10 consecutive windows
//	int const minInversionLength = 20;
//
//	//Half of window size used to compute NM in SWCPU.cpp
//	int const inversionPositionShift = 0;
//
//	int const intervalCheckLength = 50;
//
//	//Peaks (isInversion(nm) == true) that are less than maxDistance apart will be merged
//	int const maxDistance = 1;
//	int distance = maxDistance;
//
//	int startInv = -1;
//	int stopInv = -1;
//
//	int startInvRead = -1;
//	int stopInvRead = -1;
//
//	int alignStartOnRef = align.firstPosition.refPosition;
//	int alignStartOnRead = align.firstPosition.readPosition;
//
//	int alignStopOnRef = align.lastPosition.refPosition;
//	int alignStopOnRead = align.lastPosition.readPosition;
//
////TODO: improve
////	int treshold = result.NM * 4;
////	int treshold = 3;
//
////	fprintf(stderr, "Threshold: %d\n", treshold);
//
//	if (stdoutErrorProfile) {
//		int len = 0;
//		for (int i = 0; i < align.alignmentLength; ++i) {
//			printf("%s\t%llu\t%d\t%s\n", SequenceProvider.GetRefName(seqLoc.getrefId(), len), seqLoc.m_Location + align.nmPerPosition[i].refPosition, align.nmPerPosition[i].nm, readName);
//		}
//	}
//
//	return false;
//	// Only for debugging!
//	int inversionNumber = 0;
//	int len = 0;
//	for (int i = 0; i < align.alignmentLength; ++i) {
//		float nm = (32 - align.nmPerPosition[i].nm) / 32.0f;
//
//		if (nm > 0) {
////			fprintf(stderr, "%d: %d (%d)\n", positionsInRead[i], nm, treshold);
////			printf("%s\t%llu\t%llu\t%d\n",
////					SequenceProvider.GetRefName(seqLoc.getrefId(), len),
////					seqLoc.m_Location + i, seqLoc.m_Location + i + 1, nm);
//		}
//		if (startInv == -1) {
//			if (isInversion(nm)) {
//				startInv = align.nmPerPosition[i].refPosition
//						- inversionPositionShift;
//				startInvRead = align.nmPerPosition[i].readPosition
//						- inversionPositionShift;
//				stopInv = align.nmPerPosition[i].refPosition
//						- inversionPositionShift;
//				stopInvRead = align.nmPerPosition[i].readPosition
//						- inversionPositionShift;
//			}
//
//		} else {	// if(stopInv == -1) {
//			if (isInversion(nm)) {
//				stopInv = align.nmPerPosition[i].refPosition
//						- inversionPositionShift;
//				stopInvRead = align.nmPerPosition[i].readPosition
//						- inversionPositionShift;
//			} else {
//				//if (nm > 0) {
//				//					startInv = -1;
//				if (distance == 0) {
//
//					if (pacbioDebug) {
//						fprintf(stderr,
//								"Inversion detected: %d - %d, %d - %d (length: %d)\n",
//								startInv, stopInv, startInvRead, stopInvRead,
//								abs(stopInv - startInv));
//					}
//					if (abs(stopInv - startInv) > minInversionLength) {
//						if (pacbioDebug) {
//							Log.Message("Check inversion on read %s", readName);
//						}
//						inversionNumber += 1;
//						if (pacbioDebug) {
//							//Print BED files with potential inversion
//							//printf("%s\t%llu\t%llu\n", SequenceProvider.GetRefName(seqLoc.getrefId(), len), seqLoc.m_Location + startInv, seqLoc.m_Location + stopInv + 1);
//							fprintf(stderr, "Length: %d\n",
//									abs(stopInv - startInv));
//						}
//						//Positions in read and ref midpoint of inversion
//						uloc ref = (startInv + stopInv) / 2 + intervalCheckLength / 2;
//						//uloc ref = startInv;
//
//						uloc read = (startInvRead + stopInvRead) / 2;
////						uloc read = startInvRead;
//
//						char * refSeq = new char[250];
//						char * readSeq = new char[600];
//						char * revreadSeq = new char[600];
//						memset(refSeq, '\0', 250 * sizeof(char));
//						memset(readSeq, '\0', 500 * sizeof(char));
//						memset(revreadSeq, '\0', 500 * sizeof(char));
//
//						SequenceLocation testLoc;
////						testLoc.m_Location = start + align.PositionOffset + ref
////						- intervalCheckLength / 2;
//						testLoc.m_Location = start + align.PositionOffset + ref;
//						// Extract sequence from reference
//						SequenceProvider.DecodeRefSequence(refSeq, 0, testLoc.m_Location, intervalCheckLength);
//
//						// Extract sequence from read
////						strncpy(readSeq, fullReadSeq + read - intervalCheckLength, intervalCheckLength * 2);
//						strncpy(readSeq, fullReadSeq + read, intervalCheckLength * 2);
//						// Reverse complement
//						computeReverseSeq(readSeq, revreadSeq, intervalCheckLength * 2);
//
//						if (pacbioDebug) {
//							Log.Message("Ref:     %s", refSeq);
//							Log.Message("Read:    %s", readSeq);
//							Log.Message("RevRead: %s", revreadSeq);
//						}
//						if(printInvCandidateFa) {
//							printf(">%s_%d/1\n%s\n", readName, inversionNumber, refSeq);
//							printf(">%s_%d/2\n%s\n", readName, inversionNumber, revreadSeq);
//
//						}
//
//						IAlignment * aligner = new EndToEndAffine();
//
//						float scoreFwd = 0.0f;
//						float scoreRev = 0.0f;
//
//						static const float minScore = Config.GetFloat(
//								MATCH_BONUS) * 1.0f * intervalCheckLength / 4.0f;
//						aligner->SingleScore(10, 0, readSeq, refSeq, scoreFwd,
//								0);
//						aligner->SingleScore(10, 0, revreadSeq, refSeq,
//								scoreRev, 0);
//
//						delete[] refSeq;
//						refSeq = 0;
//						delete[] readSeq;
//						readSeq = 0;
//						delete[] revreadSeq;
//						revreadSeq = 0;
//						delete aligner;
//						aligner = 0;
//
//						if (pacbioDebug) {
//							Log.Message("Fwd: %f, Rev: %f, min: %f", scoreFwd, scoreRev, minScore);
//							SequenceProvider.convert(testLoc);
//						}
//
//						if(stdoutInversionBed) {
//							printf("%s\t%llu\t%llu\t%s\t%d\n", SequenceProvider.GetRefName(testLoc.getrefId(), len), testLoc.m_Location, testLoc.m_Location + intervalCheckLength, readName, 0);
//						}
//						if (scoreRev > scoreFwd && scoreRev > minScore) {
//							//if (scoreRev > scoreFwd || scoreFwd < minScore) {
//
//							inversionFound = true;
//							if (pacbioDebug) {
//								Log.Message("Inversion detected: %llu, %llu", ref,
//										read);
//								Log.Message("READ LENGTH: %d", readLength);
//								//				printf("%s\t%llu\t%llu\t%s\t%d\n", SequenceProvider.GetRefName(testLoc.getrefId(), len), testLoc.m_Location, testLoc.m_Location + 200, readName, 100);
//							}
//
//							if (interval.isReverse) {
//
//								leftOfInv.onReadStop = readLength
//								- (interval.onReadStart
//										+ alignStartOnRead);
//
//								leftOfInv.onRefStart = start
//								+ align.PositionOffset
//								+ alignStartOnRef;
//
//								leftOfInv.onReadStart = readLength
//								- (interval.onReadStart + read);
//								leftOfInv.onRefStop = start
//								+ align.PositionOffset + ref;
//
//								leftOfInv.isReverse = interval.isReverse;
//
//								rightOfInv.onReadStop = readLength
//								- (interval.onReadStart + read);
//								rightOfInv.onRefStart = start
//								+ align.PositionOffset + ref;
//								rightOfInv.onReadStart = readLength
//								- (interval.onReadStart
//										+ alignStopOnRead);
//								rightOfInv.onRefStop = start
//								+ align.PositionOffset + alignStopOnRef;
//
//								rightOfInv.isReverse = interval.isReverse;
//
//								inv.onReadStart = leftOfInv.onReadStart;
//								inv.onReadStop = leftOfInv.onReadStart;
//								inv.onRefStart = leftOfInv.onRefStop;
//								inv.onRefStop = leftOfInv.onRefStop;
//								inv.isReverse = !rightOfInv.isReverse;
//							} else {
//
//								leftOfInv.onReadStart = interval.onReadStart
//								+ alignStartOnRead;
//								leftOfInv.onRefStart = start
//								+ align.PositionOffset
//								+ alignStartOnRef;
//								leftOfInv.isReverse = interval.isReverse;
//								leftOfInv.onReadStop = interval.onReadStart
//								+ read;
//								leftOfInv.onRefStop = start
//								+ align.PositionOffset + ref;
//
//								rightOfInv.onReadStart = interval.onReadStart
//								+ read;
//								rightOfInv.onRefStart = start
//								+ align.PositionOffset + ref;
//								rightOfInv.isReverse = interval.isReverse;
//								rightOfInv.onReadStop = interval.onReadStart
//								+ alignStopOnRead;
//								rightOfInv.onRefStop = start
//								+ align.PositionOffset + alignStopOnRef;
//
//								inv.onReadStart = leftOfInv.onReadStop;
//								inv.onReadStop = leftOfInv.onReadStop;
//								inv.onRefStart = leftOfInv.onRefStop;
//								inv.onRefStop = leftOfInv.onRefStop;
//								inv.isReverse = !rightOfInv.isReverse;
//							}
//
////							if (pacbioDebug) {
////								SequenceLocation seqLoc2;
////								seqLoc2.m_Location = leftOfInv.onRefStart;
////								SequenceProvider.convert(seqLoc2);
////								Log.Message("Mapping Pos2: %llu - %llu (alignStartOnREf: %llu, corridor: %d)", leftOfInv.onRefStart, seqLoc2.m_Location, alignStartOnRef, corridor);
////								SequenceLocation seqLoc5;
////								seqLoc5.m_Location = rightOfInv.onRefStart;
////								SequenceProvider.convert(seqLoc5);
////								Log.Message("Mapping Pos5: %llu - %llu", rightOfInv.onRefStart, seqLoc5.m_Location);
////								SequenceLocation seqLoc3;
////								seqLoc3.m_Location = rightOfInv.onRefStart;
////								SequenceProvider.convert(seqLoc3);
////								Log.Message("Mapping Pos3: %llu - %llu", rightOfInv.onRefStart, seqLoc3.m_Location);
////								SequenceLocation seqLoc4;
////								seqLoc4.m_Location = rightOfInv.onRefStop;
////								SequenceProvider.convert(seqLoc4);
////								Log.Message("Mapping Pos4: %llu - %llu", rightOfInv.onRefStop, seqLoc4.m_Location);
////							}
//							return true;
//						} else {
//							//			printf("%s\t%llu\t%llu\t%s\t%d\n", SequenceProvider.GetRefName(testLoc.getrefId(), len), testLoc.m_Location, testLoc.m_Location + 200, readName, 0);
//						}
//					}
//					startInv = -1;
//					stopInv = -1;
//					distance = maxDistance;
//				} else {
//					distance -= 1;
//				}
////					getchar();
//
//			}
//		}
//	}
//
////}
////
////	extData[edIndex++] = -1;
////
////	delete[] positionsInRead;
////	positionsInRead = 0;
////
////	delete[] nmPerPos;
////	nmPerPos = 0;
//
////	getchar();
//
////	return inversionFound;
//	return false;
//}
