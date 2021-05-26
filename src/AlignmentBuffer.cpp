/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <float.h> //Eclipse
#include <limits.h>

#include "Timing.h"
#include "SequenceProvider.h"
#include "StrippedSW.h"
#include "OutputReadBuffer.h"

//#include "MemCheck.h"

using std::vector;

bool AlignmentBuffer::first = true;


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

CorridorLine * getCorridorOriginal(int const corridor, char const * readSeq,
		int & corridorHeight) {
	int const corridorWidth = corridor;
	int const qryLen = strlen(readSeq);

	corridorHeight = qryLen;
	CorridorLine * corridorLines = new CorridorLine[corridorHeight];

	NGM.Stats->corridorLen += corridorWidth;
	for (int i = 0; i < corridorHeight; ++i) {
		corridorLines[i].offset = i; // - corridorWidth / 2;
		corridorLines[i].length = corridorWidth;
	}
	return corridorLines;
}

CorridorLine * getCorridorLinear(int const corridor, char const * readSeq,
		int & corridorHeight) {
	int const corridorWidth = corridor;
	int const qryLen = strlen(readSeq);

	corridorHeight = qryLen;
	CorridorLine * corridorLines = new CorridorLine[corridorHeight];

	NGM.Stats->corridorLen += corridorWidth;
	for (int i = 0; i < corridorHeight; ++i) {
		corridorLines[i].offset = i - corridorWidth / 2;
		corridorLines[i].length = corridorWidth;
	}
	return corridorLines;
}

CorridorLine * getCorridorFull(int const corridor, char const * readSeq,
		int & corridorHeight) {
	int const corridorWidth = corridor;
	int const qryLen = strlen(readSeq);

	corridorHeight = qryLen;
	CorridorLine * corridorLines = new CorridorLine[corridorHeight];

	NGM.Stats->corridorLen += corridorWidth;
	for (int i = 0; i < corridorHeight; ++i) {
		/**
		 * 20% corridor width added on both sides for
		 * ConvexAlign "out of corridor" detection
		 * (Checks if backtracking path goes into
		 * the left/right most 10% of the corridor)
		 * TODO: fix in ConvexAlign.cpp
		 */
		corridorLines[i].offset = (int)(corridorWidth * -0.2);
		corridorLines[i].length = corridorWidth + (int)(corridorWidth * 0.2);
	}
	return corridorLines;
}

CorridorLine * getCorridorEndpoints(Interval const * interval,
		int const corridor, char const * refSeq, char const * readSeq,
		int & corridorHeight, bool const realign) {
	int corridorWidth = corridor / (realign ? 1 : 4);
	int const qryLen = strlen(readSeq);
	int const refLen = strlen(refSeq);

	corridorHeight = qryLen;
	CorridorLine * corridorLines = new CorridorLine[corridorHeight];

	float k = qryLen * 1.0f / refLen;
	float d = corridorWidth / 2.0f;

	NGM.Stats->corridorLen += corridorWidth;
	for (int i = 0; i < corridorHeight; ++i) {
		corridorLines[i].offset = (i - d) / k;
		corridorLines[i].length = corridorWidth;
	}

	return corridorLines;
}

CorridorLine * AlignmentBuffer::getCorridorEndpointsWithAnchors(Interval const * interval,
		int const corridorMultiplier, char const * refSeq, char const * readSeq,
		int & corridorHeight, int const externalQStart,
		int const readPartLength, int const fullReadLength,
		bool const realign) {
	int const qryLen = strlen(readSeq);
	int const refLen = strlen(refSeq);

	/*
	 * Find optimal corridor
	 */
	float corridorLeft = 0.0f;
	float corridorRight = 0.0f;

	float k_align = qryLen * 1.0f / refLen;
	float d_align = 0;

	for (int i = 0; i < interval->anchorLength; ++i) {

		int anchor_x = 0;
		int anchor_y = 0;
		if (interval->anchors[i].isReverse) {
			anchor_x = interval->anchors[i].onRef - interval->onRefStart;
			anchor_y = fullReadLength - interval->anchors[i].onRead
					- (readPartLength) - externalQStart;
		} else {
			anchor_x = interval->anchors[i].onRef - interval->onRefStart;
			anchor_y = interval->anchors[i].onRead - externalQStart;
		}

//		Log.Message("Anchor %d: %d, %d", i, anchor_x, anchor_y);

		float x_found = anchor_x;
		float x_expect = (anchor_y - d_align) / k_align;

//		Log.Message("x_found: %f - x_expect: %f", x_found, x_expect);

		float diff = x_expect - x_found;
		if (diff > 0) {
			corridorRight = std::max(corridorRight, diff);
		} else {
			corridorLeft = std::max(corridorLeft, diff * -1.0f);
		}
//		Log.Message("%f  -  %f", corridorLeft, corridorRight);
	}

	corridorHeight = qryLen;
	CorridorLine * corridorLines = new CorridorLine[corridorHeight];

	corridorLeft += 128;
	corridorRight += 128;

	corridorLeft = corridorLeft + (corridorLeft + corridorRight) * 0.1f;
	corridorRight = corridorRight + (corridorLeft + corridorRight) * 0.1f;

	corridorLeft = corridorLeft * corridorMultiplier;
	corridorRight = corridorRight * corridorMultiplier;

	int corridorWidth = corridorLeft + corridorRight;
	NGM.Stats->corridorLen += corridorWidth;
	for (int i = 0; i < corridorHeight; ++i) {
		corridorLines[i].offset = ((i - d_align) / k_align) - corridorRight;
		corridorLines[i].length = corridorWidth;
	}

	verbose(0, true, "Corridor length estimated from anchors is: %d", corridorWidth);

	return corridorLines;
}

char const * const AlignmentBuffer::extractReferenceSequenceForAlignment(Interval const*& interval, int & refSeqLength) {
	return extractReferenceSequenceForAlignment(interval->onRefStart, interval->onRefStop, refSeqLength);
}

char const * const AlignmentBuffer::extractReferenceSequenceForAlignment(loc const onRefStart, loc const onRefStop, int & refSeqLength) {
	char * refSeq = 0;
	if(onRefStart >= onRefStop) {
		Log.Message("Could not decode reference for alignment: %llu >= %llu.", onRefStart, onRefStop);
		return 0;
	}
	refSeqLength = onRefStop - onRefStart + 1;
	if (refSeqLength > 0) {
		//TODO: check why decoded/SequenceProvider writes outside of refSeqLen: This makes + 100 necessary
		refSeq = new char[refSeqLength + 100];
		//decode reference refSeqLength
		verbose(0, true, "RefStart: %llu, RefLength: %d", onRefStart, refSeqLength);
		if (!SequenceProvider.DecodeRefSequenceExact(refSeq, onRefStart, refSeqLength, 0)) {
			//Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
			Log.Warning("Could not decode reference for alignment");
			delete[] refSeq;
			refSeq = 0;
		}
	}
	return refSeq;
}


Align * AlignmentBuffer::computeAlignment(Interval const * interval,
		int corridor, char const * const readSeq, size_t const readLength,
		int const externalQStart, int const externalQEnd, int fullReadLength,
		MappedRead const * const read, bool realign, bool const fullAlignment,
		bool const shortRead = false) {

	if (pacbioDebug) {
		Log.Message("Read name: %s", read->name);
		Log.Message("Alignment type: realign %d, full %d, shortRead %d", realign, fullAlignment, shortRead);
	}

	if(readSeq == nullptr) {
		return 0;
	}

	static int alignmentId = 0;

	bool validAlignment = false;

	Align * align = new Align();
	try {

	#ifdef TEST_ALIGNER
		Align * alignFast = new Align();
	#endif

		int refSeqLen = 0;
		char const * refSeq = extractReferenceSequenceForAlignment(interval, refSeqLen);

		int alignedBp = 0;

		if (refSeq != 0) {

			int retryCount = 5;
			if (fullAlignment) {
				retryCount = 1;
			}

			// Corridor > refSeqLen * 2 doesn't make sense -> full matrix is computed already
			int const maxCorridorSize = refSeqLen * 2;
			corridor = std::min(corridor, maxCorridorSize);

			int const minAlignedBp = 0.05f * std::min((int) readLength, refSeqLen);

			//initialize arrays for CIGAR and MD string
			align->maxBufferLength = readLength * 4;
			align->maxMdBufferLength = readLength * 4;
			align->pBuffer1 = new char[align->maxBufferLength];
			align->pBuffer2 = new char[align->maxMdBufferLength];
			align->pBuffer1[0] = '\0';
			align->pBuffer2[0] = '\0';
			align->nmPerPostionLength = (readLength + 1) * 2;
			align->nmPerPosition = new PositionNM[align->nmPerPostionLength];

	#ifdef TEST_ALIGNER
			alignFast->pBuffer1 = new char[readLength * 4];
			alignFast->pBuffer2 = new char[readLength * 4];
			alignFast->pBuffer1[0] = '\0';
			alignFast->pBuffer2[0] = '\0';
			alignFast->nmPerPostionLength = (readLength + 1) * 2;
			alignFast->nmPerPosition = new PositionNM[alignFast->nmPerPostionLength];

	#endif


			int corridorMultiplier = 1;
			while (!validAlignment
					&& (corridor * corridorMultiplier) <= maxCorridorSize
					&& retryCount-- > 0) {

				//Local alignment
				if (pacbioDebug) {
					Log.Message("Aligning %d bp to %d bp", readLength, refSeqLen);
					Log.Message("Ref: %.*s ... %.*s", 250, refSeq, 250, refSeq + refSeqLen - 250);
					Log.Message("Read: %.*s ... %.*s", 250, readSeq, 250, readSeq - readLength - 250);
				}

				int corridorHeight = 0;
				//			CorridorLine * corridorLines = getCorridorLinear(corridor, readSeq,
				//					corridorHeight);
				// When realigning currently there are no anchors available -> larger corridor
				//TODO: Use anchors from original interval, or extract anchors from first alignment!
				CorridorLine * corridorLines = 0;
				if(fullAlignment) {
					corridorLines = getCorridorFull(refSeqLen, readSeq,
							corridorHeight);
				} else {
					if(shortRead) {
						verbose(0, true, "Corridor width: %d", corridor * corridorMultiplier);
						corridorLines = getCorridorLinear(corridor * corridorMultiplier, readSeq,
								corridorHeight);
					} else {
						if(corridorMultiplier < 3 && !realign && interval->anchorLength > 0) {
							corridorLines = getCorridorEndpointsWithAnchors(interval,
									corridorMultiplier, refSeq, readSeq, corridorHeight, externalQStart, readPartLength, fullReadLength, realign);
						} else {
							verbose(0, true, "Corridor width: %d", corridor * corridorMultiplier);
							corridorLines = getCorridorEndpoints(interval,
									corridor * corridorMultiplier, refSeq, readSeq, corridorHeight, realign);
						}
					}
				}

				Timer algnTimer;
				algnTimer.ST();

				if (stdoutPrintAlignCorridor) {

					for(int x = 0; x < interval->anchorLength; ++x) {
						if(interval->anchors[x].isReverse) {
							printf("%d\t%d\t%lld\t%d\t%d\n", alignmentId, read->ReadId, interval->anchors[x].onRef - interval->onRefStart,
									fullReadLength - interval->anchors[x].onRead - (readPartLength) - externalQStart, 3);
						} else {
							printf("%d\t%d\t%lld\t%d\t%d\n", alignmentId, read->ReadId, interval->anchors[x].onRef - interval->onRefStart,
									interval->anchors[x].onRead - externalQStart, 3);
						}
					}
					printf("%d\t%d\t%d\t%s\t%d\n", alignmentId, read->ReadId, read->ReadId,
							read->name, -4);
					printf("%d\t%d\t%d\t%d\t%d\n", alignmentId, read->ReadId, interval->isReverse,
							corridorLines[0].length, -5);
					printf("%d\t%d\t%d\t%d\t%d\n", alignmentId, read->ReadId, externalQStart,
							externalQEnd, -6);

				}

				if(pacbioDebug) {
					int cellCount = 0;
					for(int i = 0; i < corridorHeight; ++i) {
						cellCount += corridorLines[i].length;
					}
					Log.Message("Computing %d cells of alignment Matrix", cellCount);
				}

				//Hack to pass readId to convex alignment class for plotting
				//TODO: remove
				align->svType = read->ReadId;

				#ifdef TEST_ALIGNER
					alignFast->svType = read->ReadId;
				#endif

	#ifdef TEST_ALIGNER
				Timer tmr1;
				tmr1.ST();
	#endif


				if(pacbioDebug) {
					Log.Message("ExternalQstart: %d, ExternalQEnd: %d", externalQStart, externalQEnd);
				}
				int cigarLength = aligner->SingleAlign(alignmentId, corridorLines,
						corridorHeight, (char const * const ) refSeq,
						(char const * const ) readSeq, *align, externalQStart, externalQEnd, 0);
//				cigarLength = -1;

	#ifdef TEST_ALIGNER
				float time1 = tmr1.ET();
				Timer tmr2;
				tmr2.ST();
				int const cigarLengthFast = alignerFast->SingleAlign(alignmentId, corridorLines,
						corridorHeight, (char const * const ) refSeq,
						(char const * const ) readSeq, *alignFast, externalQStart, externalQEnd, 0);

				float time2 = tmr2.ET();

				Log.Message("%d/%d bp: %f - %f", strlen(refSeq), strlen(readSeq), time1, time2);

				if(!(cigarLengthFast == cigarLength && alignFast->Score == align->Score && strcmp(align->pBuffer1, alignFast->pBuffer1) == 0 && strcmp(align->pBuffer2, alignFast->pBuffer2) == 0)) {
					Log.Message("Ref:  %s", refSeq);
					Log.Message("Read: %s", readSeq);
					Log.Message("Convex:     %d %f %s", cigarLength, align->Score, align->pBuffer1);
					Log.Message("ConvexFast: %d %f %s", cigarLengthFast, alignFast->Score, alignFast->pBuffer1);
					Log.Error("Not equal");
				}
	#endif

				alignmentId += 1;

				if (pacbioDebug) {
					Log.Message("Aligning took %f seconds", algnTimer.ET());
					Log.Message("CIGAR: %.*s", 250, align->pBuffer1);
					Log.Message("MD:    %.*s", 250, align->pBuffer2);
				}
				delete[] corridorLines;
				corridorLines = 0;

				alignedBp = readLength - (align->QStart - externalQStart) - (align->QEnd - externalQEnd);

				validAlignment = cigarLength == fullReadLength; // && alignedBp > minAlignedBp;

				if(!validAlignment) {
					//corridor = corridor * 2;
					corridorMultiplier += 1;
					if (pacbioDebug) {
						Log.Message("Invalid alignment found. Running again with corridor %d, %d attempts left", corridor * corridorMultiplier, retryCount);
					}
					NGM.Stats->invalidAligmentCount += 1;
				}
			}
			delete[] refSeq;
			refSeq = 0;

			#ifdef TEST_ALIGNER
				delete alignFast; alignFast = 0;
			#endif
		} else {
			Log.Error("Could not extract reference sequence for read %s.", read->name);
			validAlignment = false;
		}

		if (validAlignment) {
			if (pacbioDebug) {
				Log.Message("%d of %d bp successfully aligned with score %f and identity %f", alignedBp, readLength, align->Score, align->Identity);
			}
			NGM.Stats->alignmentCount += 1;
		} else {
			if (pacbioDebug) {
				Log.Message("Could not align sequences.");
			}
			// If alignment failed delete align object and return 0
			if(align != 0) {
				align->clearBuffer();
				align->clearNmPerPosition();
				delete align;
				align = 0;
			}
		}
	} catch (...) {
		Log.Message("Warning: could not compute alignment for read %s", read->name);
		// If alignment failed delete align object and return 0
		if(align != 0) {
			align->clearBuffer();
			align->clearNmPerPosition();
			delete align;
			align = 0;
		}
	}
	return align;
}

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

				loc iRef = anchors[i].onRef;
				loc jRef = anchors[j].onRef;

				loc refDiff = 0;
				if (anchors[j].isReverse) {
					refDiff = jRef - iRef;
				} else {
					refDiff = iRef - jRef;
				}

				loc readDiff = anchors[i].onRead - anchors[j].onRead;

				loc diff = llabs(refDiff - readDiff);
				loc maxDiff = std::max(llabs(refDiff), readDiff) * 0.25;
				loc maxRefDiff = readPartLength * 2.0f;
//				Log.Message("cLIS: (rev %d, %d) ref %lld - %lld = %lld; read %d - %d = %d; diff %lld < %lld", anchors[j].isReverse, anchors[i].isReverse, iRef, jRef, refDiff, anchors[i].onRead, anchors[j].onRead, readDiff, diff, maxDiff);

				if (DP[j] + 1 > DP[i] // New LIS must be longer then best so far
				&& anchors[j].isReverse == anchors[i].isReverse //Anchors have to be on the same strand to be joined!
						&& (diff < maxDiff // The distance between the two anchors on the read and on the reference should not exceed a threshold. Adhoc threshold of 15% (INDEL sequencing error was to littel). 25% works better.
								|| (anchors[i].onRead == anchors[j].onRead
										&& llabs(refDiff) <= readPartLength) // if anchors stem from the same read position accept ref positions that are in close proximity (to avoid additional segments caused be anchors that are only a few bp off)
						) && refDiff < maxRefDiff // Distance between anchors must be smaller than twice the anchor length. With this, gaps break the segment even if they are on the same "diagonal"
						&& refDiff >= 0 // And is not allowed to be negative (only if read position is the same, see above)
								) {
//					Log.Message("---> OK!");
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

bool AlignmentBuffer::isSameDirection(Interval const * a, Interval const * b) {
	return a->isReverse == b->isReverse;
}

bool isAnchorInCorridor(REAL k, REAL d, REAL corridor, Anchor const & testee) {

	loc onRefStart = testee.onRef;

	bool inCorridor = false;
	REAL y = testee.onRead;
	loc upperRefStart = (loc) round((y - (d + corridor)) / k);
	loc lowerRefStart = (loc) round((y - (d - corridor)) / k);

	if (upperRefStart < lowerRefStart) {
		REAL tmp = upperRefStart;
		upperRefStart = lowerRefStart;
		lowerRefStart = tmp;
	}

	inCorridor = inCorridor
			|| (onRefStart <= upperRefStart && onRefStart >= lowerRefStart);
//	Log.Message("Start: %llu <= %llu <= %llu", testee.onRefStart, lowerRefStart, upperRefStart);

	return inCorridor;
}

bool AlignmentBuffer::isIntervalInCorridor(REAL k, REAL d, REAL corridor,
		Interval const * testee, bool const switched) {

	loc onRefStart = testee->onRefStart;
	loc onRefStop = testee->onRefStop;

	if (switched) {
		onRefStart = testee->onRefStop;
		onRefStop = testee->onRefStart;
	}

	bool inCorridor = false;
	REAL y = testee->onReadStart;
	loc upperRefStart = (loc) round((y - (d + corridor)) / k);
	loc lowerRefStart = (loc) round((y - (d - corridor)) / k);

	if (upperRefStart < lowerRefStart) {
		REAL tmp = upperRefStart;
		upperRefStart = lowerRefStart;
		lowerRefStart = tmp;
	}

	verbose(3, true, "%llu <= %llu && %llu >= %llu", onRefStart, upperRefStart, onRefStart, lowerRefStart);
	verbose(3, true, "onRefStart <= upperRefStart = %d, onRefStart >= lowerRefStart = %d", onRefStart <= upperRefStart, onRefStart >= lowerRefStart);
	inCorridor = inCorridor
			|| (onRefStart <= upperRefStart && onRefStart >= lowerRefStart);
//	Log.Message("Start: %llu <= %llu <= %llu", testee.onRefStart, lowerRefStart, upperRefStart);

	y = testee->onReadStop;
	loc upperRefStop = (loc) round((y - (d + corridor)) / k);
	loc lowerRefStop = (loc) round((y - (d - corridor)) / k);

	if (upperRefStop < lowerRefStop) {
		REAL tmp = upperRefStop;
		upperRefStop = lowerRefStop;
		lowerRefStop = tmp;
	}

	inCorridor = inCorridor
			&& (onRefStop <= upperRefStop && onRefStop >= lowerRefStop);
	verbose(3, true, "%llu <= %llu && %llu >= %llu", onRefStop, upperRefStop, onRefStop, lowerRefStop);
	verbose(3, true, "onRefStop <= upperRefStop = %d, onRefStop >= lowerRefStop = %d", onRefStop <= upperRefStop, onRefStop >= lowerRefStop);
//	Log.Message("Stop: %llu <= %llu <= %llu", testee.onRefStop, lowerRefStop, upperRefStop);

//TODO: add max distance on read

	return inCorridor;
}

Interval * AlignmentBuffer::toInterval(Anchor const & anchor) {
//	if (intervalBufferIndex >= 1000) {
//		return 0;
//	}

	Interval * interval = new Interval();

//	intervalBuffer[intervalBufferIndex++] = interval;

	interval->isReverse = anchor.isReverse;
	interval->score = anchor.score;
	interval->onReadStart = interval->onReadStop = anchor.onRead;
	interval->onRefStart = interval->onRefStop = anchor.onRef;

	interval->anchorLength = 1;
	interval->anchors = new Anchor[interval->anchorLength];
	interval->anchors[0] = anchor;

	// To avoid that single anchor intervals that didn't pass cLIS
	// start a new mapped segment
	interval->isAssigned = true;

	return interval;
}

void AlignmentBuffer::addAnchorAsInterval(Anchor const & anchor,
		MappedSegment & segment) {
	if (segment.length < 1000) {
		Interval * interval = toInterval(anchor);
		if (interval != 0) {
			segment.list[segment.length++] = interval;
		}
	}
}

//bool AlignmentBuffer::isCompatible(Anchor const & anchor,
//		Interval const * interval) {
//	bool isCompatible = false;
//
//	REAL corridorSize = 512;
//
//	if (anchor.isReverse == interval->isReverse) {
////		isCompatible = (anchor.onRead < interval->onReadStart
////				|| anchor.onRead > interval->onReadStop)
////				&& (anchor.onRef < interval->onRefStart
////						|| anchor.onRef > interval->onRefStop);
//
//		isCompatible = isAnchorInCorridor(interval->m, interval->b,
//				corridorSize, anchor);
//	}
//
//	return isCompatible;
//}

//bool AlignmentBuffer::isCompatible(Anchor const & anchor,
//		MappedSegment const & segment) {
//	bool compatible = false;
//
//	for (int i = 0; i < segment.length && !compatible; ++i) {
//
//		compatible = isCompatible(anchor, segment.list[i]);
//	}
//
//	return compatible;
//}

// Checks whether a is compatible (is located in the "corridor"
// of b.
bool AlignmentBuffer::isCompatible(Interval const * a, Interval const * b, REAL corridorSize) {

	bool isCompatible = false;
	// Trade off (at the moment):
	// Bigger corridor: higher chance that alignment won't span event (score too low)
	// and a part of the reads is not aligned
	// Smaller corridor: reads are split, SVs are not visible in one alignment
	// TODO: after alignment add check if whole read was aligned. if not,
	// realign unanligned part
//	REAL corridorSize = 2048 * 4;

	// Check if regression worked on read
	if (b->m != 0 && b->b != 0 && (b->r * b->r) > 0.8f) {

		if (a->isReverse == b->isReverse) {
			// a and b are on the same strand
			isCompatible = isIntervalInCorridor(b->m, b->b, corridorSize, a,
					false);
		} else {
			if (pacbioDebug) {
				Log.Message("Intervals not on same strand");
			}
			/*
			 * a and b are on different strands
			 * Check if reads spans an inversion
			 * If read spans an inversion, the
			 * inverted part must be added to the segment
			 * in order to keep consolidateSegements
			 * to join the non inverted parts
			 * through the inversion.
			 */

			// Switches the direction of the testee
			isCompatible = isIntervalInCorridor(b->m, b->b, corridorSize, a, true);

			// Switches the direction of the tester and uses the testeee regression
			// to check whether both are compatible
			isCompatible = isCompatible || isIntervalInCorridor(a->m, a->b, corridorSize, b, true);
		}

	}

	return isCompatible;
}

bool AlignmentBuffer::canSpanDeletionInsertion(Interval const * a, Interval const * b, REAL corridorSize) {

	bool merge = true;

	int distanceOnRead = getDistanceOnRead(a, b);
	loc distanceOnRef = getDistanceOnRef(a, b);

	verbose(3, true, "DistOnRef: %llu, DistOnRead: %d, Corridor: %f", distanceOnRef, distanceOnRead, corridorSize);
	if(distanceOnRead < (2 * readPartLength)) {
		//Possible deletion
		merge = abs(distanceOnRef - distanceOnRead) < corridorSize;
	} else if(distanceOnRef < (2 * readPartLength)) {
		//Possible insertion
		merge = abs(distanceOnRead - distanceOnRef) < corridorSize;
	} else {
		/**
		 * Do nothing, intervals show a larger gap. Probably, poor quality read.
		 * Better to merge than to risk fragmented alignments
		 */
		merge = abs(distanceOnRead - distanceOnRef) < corridorSize;
	}
	return merge;
}

bool AlignmentBuffer::spansChromosomeBorder(Interval const * a, Interval const * b) {

	_SequenceProvider::Chromosome chrA = SequenceProvider.getChrStart((a->onRefStop + a->onRefStart) / 2);
	_SequenceProvider::Chromosome chrB = SequenceProvider.getChrStart((b->onRefStop + b->onRefStart) / 2);

	verbose(0, true, "#### Intervals on same chr:");
	verbose(0, true, "#### 1: %llu", (a->onRefStop + a->onRefStart) / 2);
	verbose(0, true, "#### 2: %llu", (b->onRefStop + b->onRefStart) / 2);
	verbose(0, true, "#### Chr1: %llu - %llu", chrA.start, chrA.end);
	verbose(0, true, "#### Chr2: %llu - %llu", chrB.start, chrB.end);

	return chrA.start != chrB.start;
}

bool AlignmentBuffer::isContained(Interval const * a, Interval const * b) {
	return a->onReadStart >= b->onReadStart && a->onReadStop <= b->onReadStop

	&& a->onRefStart >= b->onRefStart && a->onRefStop <= b->onRefStop
			&& a->isReverse == b->isReverse;
}


Interval * AlignmentBuffer::mergeIntervals(Interval * a, Interval * b) {

	if (a->onReadStart > b->onReadStart) {
		a->onReadStart = b->onReadStart;
		a->onRefStart = b->onRefStart;

	}
	if (a->onReadStop < b->onReadStop) {
		a->onReadStop = b->onReadStop;
		a->onRefStop = b->onRefStop;
	}

	a->score += b->score;

	int anchorCount = a->anchorLength + b->anchorLength;
	Anchor * tmpA = a->anchors;

	a->anchors = new Anchor[anchorCount];
	memcpy(a->anchors, tmpA, a->anchorLength * (sizeof(Anchor)));
	memcpy(a->anchors + a->anchorLength, b->anchors,
			b->anchorLength * (sizeof(Anchor)));
	a->anchorLength = anchorCount;
	delete[] tmpA;
	tmpA = 0;

	a->isAssigned = a->isAssigned && b->isAssigned;

	return a;
}

/**
 * Overlap on ref coordinates >= readPartLength
 * Overlap on read coordinadtes <= readPartLength
 * Difference in overlap must be > 0
 * a and b must be on same strand!
 */
bool AlignmentBuffer::isDuplication(Interval const * a,
		Interval const * b, loc & dupLength) {

	verbose(3, true, "isDuplication:");
	verbose(3, "", a);
	verbose(3, "", b);

	int overlapOnRead = getOverlapOnRead(a, b);
	loc overlapOnRef = 0;

	if (a->isReverse) {
		// if a and b are on reverse strand, refstart and refstop must be switched
		overlapOnRef = std::max(0ll, std::min(a->onRefStart, b->onRefStart) - std::max(a->onRefStop, b->onRefStop));
	} else {
		overlapOnRef = std::max(0ll, std::min(a->onRefStop, b->onRefStop) - std::max(a->onRefStart, b->onRefStart));
	}

	verbose(2, true, "overlap on read: %d, overlap on ref: %lld", overlapOnRead, overlapOnRef);

	loc overlapDiff = std::max(0ll, overlapOnRef - overlapOnRead);
	dupLength = overlapDiff;

	return overlapOnRef >= (int) (1.0f * readPartLength)
			&& overlapOnRead <= (int) (readPartLength * 1.0f) && overlapDiff > 0ll; // overlapDiff > (int)round(readPartLength / 2.0);
}

bool sortIntervalsInSegment(Interval const * a, Interval const * b) {
	return a->onReadStart < b->onReadStart;
}

bool sortIntervalsByScore(Interval const * a, Interval const * b) {
	return a->score > b->score;
}

/**
 * - Sort sub-reads by read position
 * - Compute max. number of HSP (based on read length)
 * - Run cLIS algorithm on reference positions to retrieve HSP
 * - If isUnique, build Interval from sub-read set and compute regression
 */
Interval * * AlignmentBuffer::getIntervalsFromAnchors(int & intervalsIndex, Anchor * allFwdAnchors, int allFwdAnchorsLength, Anchor * allRevAnchors, int allRevAnchorsLength, MappedRead * read) {

	verbose(0, true, "Finding LIS for read %d", read->ReadId);

	//Sort by position on read. Probably not necessary!!
	std::sort(allFwdAnchors, allFwdAnchors + allFwdAnchorsLength,
			sortAnchorOnRead);

	// Change lo
	maxSegmentCount = std::max(10, Config.getMaxSegmentNumberPerKb(read->length) * 2);
	int const maxcLISRunNumber = maxSegmentCount;
	Interval * * intervals = new Interval * [maxcLISRunNumber];
	intervalsIndex = 0;

	int cLISRunNumber = 0;
	//int const maxRunNumber = std::min(maxcLISRunNumber, Config.getMaxCLISRuns());
	int const maxRunNumber = Config.getMaxCLISRuns();
	verbose(0, true, "maxRunNumber: %d, maxcLISRunNumber: %d, MaxCLISRuns: %d", maxRunNumber, maxcLISRunNumber, Config.getMaxCLISRuns());
	int runNumber = 0;
	bool finished = false;
	// TODO: find better stop criteria!
	while (cLISRunNumber < maxcLISRunNumber && !finished && ++runNumber < maxRunNumber) {
		if (allFwdAnchorsLength > 0) {

			//Find constrained LIS
			int lisLength = 0;
			int * lis = cLIS(allFwdAnchors, allFwdAnchorsLength, lisLength);


			if(lisLength < 1) {
				//Nothing found
				finished = true;
			} else {

				int minOnRead = INT_MAX;
				int maxOnRead = 0;

				loc minOnRef = LLONG_MAX;
				loc maxOnRef = 0;

				bool isReverse = false;

				float intervalScore = 0.0f;

				verbose(1, false, "cLIS %d:", cLISRunNumber);

				//Remove LIS from candidates
				int posInLIS = lisLength - 1;
				allRevAnchorsLength = 0;

				std::unique_ptr<REAL[]> regX(new REAL[std::max(2, allFwdAnchorsLength)]);
				std::unique_ptr<REAL[]> regY(new REAL[std::max(2, allFwdAnchorsLength)]);
				int pointNumber = 0;

				//Find intervall that covers all anchors best
				Interval * interval = new Interval();

				interval->anchors = new Anchor[lisLength + 1];
				interval->anchorLength = 0;

				/**
				 * Only intervals that have at least one
				 * unique (only on mapping position) anchor
				 * are considered.
				 */
//				bool isUnique = lisLength >= 4;
				bool isUnique = false;
				for (int i = 0; i < allFwdAnchorsLength; ++i) {
					if (posInLIS >= 0 && i == lis[posInLIS]) {

						interval->anchors[interval->anchorLength++] = allFwdAnchors[lis[posInLIS]];

						isUnique = isUnique || allFwdAnchors[lis[posInLIS]].isUnique;

						int onRead = allFwdAnchors[lis[posInLIS]].onRead;

						//Print current LIS
						verbose(0, false, "%d, ", lis[posInLIS]);

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
						} else {

							if(onRead < minOnRead) {
								minOnRead = onRead;
								minOnRef = allFwdAnchors[lis[posInLIS]].onRef;
							}

							if((onRead + readPartLength) > maxOnRead) {
								maxOnRead = onRead + readPartLength;
								maxOnRef = allFwdAnchors[lis[posInLIS]].onRef + readPartLength;
							}
						}

						regY[pointNumber] = onRead;
						if(isReverse) {
							regX[pointNumber] = allFwdAnchors[lis[posInLIS]].onRef + readPartLength;
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
				verbose(0, true, "");

				if(isUnique) {

					//Debug print anchors used in cLIS
					for(int i = 0; i < interval->anchorLength; ++i) {
						Anchor & anchor = interval->anchors[i];
						if(anchor.isReverse) {
							printDotPlotLine(read->ReadId, read->name,
									anchor.onRead,
									anchor.onRead + readPartLength,
									anchor.onRef + readPartLength,
									anchor.onRef,
									anchor.score,
									anchor.isReverse,
									DP_TYPE_CLIS + cLISRunNumber, DP_STATUS_OK);
						} else {
							printDotPlotLine(read->ReadId, read->name,
									anchor.onRead,
									anchor.onRead + readPartLength,
									anchor.onRef,
									anchor.onRef + readPartLength,
									anchor.score,
									anchor.isReverse,
									DP_TYPE_CLIS + cLISRunNumber, DP_STATUS_OK);
						}
					}


					// Linear regression for segment
					REAL m,b,r;
					if(pointNumber == 1) {
						regX[0] = minOnRef;
						regY[0] = minOnRead;
						regX[1] = maxOnRef;
						regY[1] = maxOnRead;
						pointNumber = 2;
					}
					linreg(pointNumber,regX,regY,&m,&b,&r);
					verbose(0, true, "Regression: m=%.*f b=%.*f r=%.*f\n", DBL_DIG,m,DBL_DIG,b,DBL_DIG,r);

					interval->isReverse = isReverse;
					interval->score = intervalScore;

					interval->onReadStart = minOnRead;
					interval->onReadStop = maxOnRead;

					interval->onRefStart = minOnRef;
					interval->onRefStop = maxOnRef;

					interval->m = m;
					interval->b = b;
					interval->r = r;

					if(interval->lengthOnRead() > 0 && interval->lengthOnRef() > 0ll) {
						verbose(0, true, "New interval: ");
						verbose(1, "", interval);
						intervals[intervalsIndex++] = interval;
					} else {
						if(pacbioDebug) {
							Log.Message("Skipping interval. Too short. Deleting.");
							if(interval != 0) {
								delete interval;
								interval = 0;
							}
						}
					}
					cLISRunNumber += 1;

//					if(interval->m != 0.0) {
//						printDotPlotLine(read->ReadId, read->name,
//								interval->m, interval->b, r,
//								interval->score,
//								interval->isReverse,
//								DP_TYPE_SEQMENTS_REG + cLISRunNumber, DP_STATUS_NOCOORDS);
//
//					}

					printDotPlotLine(read->ReadId, read->name,
							interval->onReadStart,
							interval->onReadStop,
							interval->onRefStart,
							interval->onRefStop,
							interval->score,
							interval->isReverse,
							DP_TYPE_SEQMENTS + cLISRunNumber, DP_STATUS_OK);

					//TODO: check chromosome borders
				} else {
					if(interval != 0) {
						delete interval;
						interval = 0;
					}
//					if(lisLength == 1) {
//						finished = true;
//					}
				}

				//Switch first with second list
				Anchor * tmp = allFwdAnchors;
				allFwdAnchors = allRevAnchors;
				allFwdAnchorsLength = allRevAnchorsLength;
				allRevAnchors = tmp;
				allRevAnchorsLength = 0;

			}

			lisLength = 0;
			delete[] lis;

		} else {
			// If no unprocessed anchors are found anymore
			finished = true;
		}
		verbose(0, true, "cLISRunNumber %d < maxcLISRunNumber %d && !finished %d && ++runNumber %d < maxRunNumber %d", cLISRunNumber, maxcLISRunNumber, !finished, runNumber, maxRunNumber);
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

int AlignmentBuffer::checkForSV(Align const * const align, Interval const * interval, char const * const fullReadSeq, uloc inversionMidpointOnRef, uloc inversionMidpointOnRead, int const inversionLength, MappedRead * read) {

	int inversionNumber = 0;

	static bool const noLowQualSplit = !Config.getLowQualitySplit();

	int const readCheckLength = 50;
	// Length of regions that will be checked left and right to the midpoint of the
	// inversion on the reference
	int const refCheckLength = 250;


	int const minInversionLength = 10;

	int svType = SV_NONE;

	if (inversionLength > minInversionLength) {
		inversionNumber += 1;

		SequenceLocation inversionCheckLocation;
		inversionCheckLocation.m_Location = interval->onRefStart + align->PositionOffset + inversionMidpointOnRef - refCheckLength - (inversionLength / 2);

//		svType = alignmentCheckForInversion(inversionLength, refCheckLength, inversionCheckLocation, inversionMidpointOnRead, read->name, inversionNumber, fullReadSeq);

		int const fullReadSeqLength = strlen(fullReadSeq);

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

		// Extract sequence from read.
		if (pacbioDebug) {
			Log.Message("CheckLength: %d, Midpoint: %d, readSeqLength: %d", readCheckLength, inversionMidpointOnRead, readSeqLenght);
		}
		if (readCheckLength <= inversionMidpointOnRead && (inversionMidpointOnRead + readCheckLength) < fullReadSeqLength) {
			strncpy(readSeq, fullReadSeq + inversionMidpointOnRead - readCheckLength, readCheckLength * 2);
		}

		if (strlen(readSeq) > 0) {
			// Reverse complement
			computeReverseSeq(readSeq, revreadSeq, readCheckLength * 2);
			if (pacbioDebug) {
				Log.Message("Ref:      %s", refSeq);
				Log.Message("Read:     %s", readSeq);
				Log.Message("RevRead:  %s", revreadSeq);
			}
			if (printInvCandidateFa) {
				printf(">%s_%d/1\n%s\n", read->name, inversionNumber, refSeq);
				printf(">%s_%d/2\n%s\n", read->name, inversionNumber, revreadSeq);
			}
			IAlignment * lqCheckAligner = new StrippedSW();

			float scoreFwd = 0.0f;
			float scoreRev = 0.0f;
			static const float minScore = 1.0f * readCheckLength / 4.0f;
			//Formula from Moritz -> p-val > 0.05 for Match:1, mismacth: -4, gap open: -2, gap extend -1
//			const float minScore = 15.0f;
			lqCheckAligner->SingleScore(10, 0, refSeq, readSeq, scoreFwd, 0);
			lqCheckAligner->SingleScore(10, 0, refSeq, revreadSeq, scoreRev, 0);

			delete lqCheckAligner;
			lqCheckAligner = 0;

			if ((scoreRev / scoreFwd) > Config.getInvScoreRatio() && scoreRev > minScore) {
//			if (scoreRev >= scoreFwd && scoreRev > minScore) {
				svType = SV_INVERSION;
			} else if (scoreRev < minScore && scoreFwd < minScore && !noLowQualSplit) {

				svType = SV_TRANSLOCATION;
			}
			//		inversionVerified = scoreRev > scoreFwd || scoreFwd < minScore;
		} else {
			verbose(0, true, "Skipping alignment. Readpart not long enough.");
		}
		delete[] refSeq;
		refSeq = 0;
		delete[] readSeq;
		readSeq = 0;
		delete[] revreadSeq;
		revreadSeq = 0;

		if (svType != SV_NONE) {
			if (pacbioDebug) {
				if (svType == SV_INVERSION) {
					Log.Message("Inversion detected: %llu, %llu", inversionMidpointOnRef,
							inversionMidpointOnRead);
				} else if(svType == SV_TRANSLOCATION) {
					Log.Message("Translocation detected: %llu, %llu", inversionMidpointOnRef,
							inversionMidpointOnRead);
				}
			}
		}
	} else {
		if(pacbioDebug) {
			Log.Message("Inversion too short!");
		}
	}
	return svType;
}

int AlignmentBuffer::detectMisalignment(Align const * const align,
		Interval const * alignedInterval, char const * const readPartSeq,
		Interval * leftOfInv, Interval * rightOfInv, MappedRead * read) {


	//*********************//
	//Detect inversions
	//*********************//

	uloc bestInversionMidpointOnRead = 0;
	uloc bestInversionMidpointOnRef = 0;

	//Half of window size used to compute NM in SWCPU.cpp
	int const inversionPositionShift = 0;

	/**
	 * Max number of peaks allowed for read.
	 * If more peaks are found, read is considered to be
	 * of bad quality and no split is performed
	 */
	int const maxCheckCount = std::max(1, (int)((read->length / 1000.0f) / 2.0f));

//Peaks (isInversion(nm) == true) that are less than maxDistance apart will be merged
	int const maxDistance = 20;
	int distance = maxDistance;

	int startInv = -1;
	int stopInv = -1;

	int startInvRead = -1;
	int stopInvRead = -1;

	SequenceLocation mappingLocation;
	mappingLocation.m_Location = alignedInterval->onRefStart
			+ align->PositionOffset;
	SequenceProvider.convert(mappingLocation);

	if (stdoutErrorProfile) {
		int len = 0;
		for (int i = 0; i < align->alignmentLength; ++i) {
			printf("%s\t%llu\t%d\t%s\n", SequenceProvider.GetRefName(mappingLocation.getrefId(), len), mappingLocation.m_Location + align->nmPerPosition[i].refPosition, align->nmPerPosition[i].nm, read->name);
		}
	}

	/**
	 * Reconds number of detected peaks
	 */
	int checkCount = 0;

	// Extremely simple "peak-finding" algorithm
	int len = 0;
	int result = SV_NONE;
	int bestResult = SV_NONE;
	for (int i = 0; i < align->alignmentLength /*&& result == SV_NONE*/; ++i) {
		float nm = (32 - align->nmPerPosition[i].nm) / 32.0f;

		//Log.Message("Inv: %d - %d", startInvRead, stopInvRead);

		if (startInv == -1) {
			if (isInversion(nm)) {
				startInv = align->nmPerPosition[i].refPosition
						- inversionPositionShift;
				startInvRead = align->nmPerPosition[i].readPosition
						- inversionPositionShift;
				stopInv = align->nmPerPosition[i].refPosition
						- inversionPositionShift;
				stopInvRead = align->nmPerPosition[i].readPosition
						- inversionPositionShift;
			}

		} else {	// if(stopInv == -1) {
			if (isInversion(nm)) {
				stopInv = align->nmPerPosition[i].refPosition
						- inversionPositionShift;
				stopInvRead = align->nmPerPosition[i].readPosition
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
						Log.Message("Low alignment quality region detected: %d - %d, %d - %d (length: %d) on read %s",
						startInv, stopInv, startInvRead, stopInvRead,
						abs(stopInv - startInv), read->name);
					}

					uloc inversionMidpointOnRef = (startInv + stopInv) / 2;
					uloc inversionMidpointOnRead = (startInvRead + stopInvRead) / 2;
					int inversionLength = abs(stopInv - startInv);

//					char * test = new char[stopInvRead - startInvRead + 100];
//					strncpy(test, readPartSeq + startInvRead, (stopInvRead - startInvRead));
//					test[(stopInvRead - startInvRead)] = '\n';
//
//					SequenceLocation inversionCheckLocation;
//					inversionCheckLocation.m_Location = alignedInterval->onRefStart + align->PositionOffset + startInv;
//					char * refSeq = new char[inversionLength + 100];
//					SequenceProvider.DecodeRefSequence(refSeq, 0, inversionCheckLocation.m_Location, inversionLength);
//
//					Log.Message("Inv read sequence: %s", test);
//					Log.Message("Inv read sequence: %s", refSeq);

					checkCount += 1;
					result = checkForSV(align, alignedInterval, readPartSeq, inversionMidpointOnRef, inversionMidpointOnRead, inversionLength, read);

					if(bestResult == SV_NONE || result == SV_INVERSION) {
						bestResult = result;
						bestInversionMidpointOnRef = inversionMidpointOnRef;
						bestInversionMidpointOnRead = inversionMidpointOnRead;
					}

					startInv = -1;
					stopInv = -1;
					startInvRead = -1;
					stopInvRead = -1;
					distance = maxDistance;
				} else {
					distance -= 1;
				}
			}
		}
	}

	if(checkCount <= maxCheckCount) {
		if (bestResult != SV_NONE) {
			if (alignedInterval->isReverse) {

				// alignStartOnRead is not relative to the full read, but to the part (interval)
				// that was aligned in the first place (the alignemnt that is scanned for an inversion)
				// Qstart is realtive to the full read.
				// However leftOfInv must be relative to the full read as well.
				// To define the stop position (start on the fwd strand) QStart can be used
				// For the start position (stop on the fwd strand) we use the position were
				// the inversion was detected. Since this position in relative to the aligned part
				// of the read (interval) we have to add additionalQstart
				int additionalQStart = align->QStart - align->firstPosition.readPosition;
				leftOfInv->onReadStop = read->length - (align->QStart);
				leftOfInv->onReadStart = read->length - (additionalQStart + bestInversionMidpointOnRead);

				// Position on the reference works perfectly with the mapping position
				// from the intial alignment. No need to account for the fact that only
				// a part of the read was aligned
				leftOfInv->onRefStart = alignedInterval->onRefStart + align->PositionOffset + align->firstPosition.refPosition;
				leftOfInv->onRefStop = alignedInterval->onRefStart + align->PositionOffset + bestInversionMidpointOnRef;

				leftOfInv->isReverse = alignedInterval->isReverse;

				// Same as for leftOfInv
				rightOfInv->onReadStart = read->length - (align->lastPosition.readPosition + additionalQStart);
				rightOfInv->onReadStop = read->length - (bestInversionMidpointOnRead + additionalQStart);

				rightOfInv->onRefStart = alignedInterval->onRefStart + align->PositionOffset + bestInversionMidpointOnRef;
				rightOfInv->onRefStop = alignedInterval->onRefStart + align->PositionOffset + align->lastPosition.refPosition;

				rightOfInv->isReverse = alignedInterval->isReverse;
			} else {
				leftOfInv->onReadStart = alignedInterval->onReadStart + align->firstPosition.readPosition;
				leftOfInv->onReadStop = alignedInterval->onReadStart + bestInversionMidpointOnRead;

				leftOfInv->onRefStart = alignedInterval->onRefStart + align->PositionOffset + align->firstPosition.refPosition;
				leftOfInv->onRefStop = alignedInterval->onRefStart + align->PositionOffset + bestInversionMidpointOnRef;
				leftOfInv->isReverse = alignedInterval->isReverse;

				rightOfInv->onReadStart = alignedInterval->onReadStart + bestInversionMidpointOnRead;
				rightOfInv->onReadStop = alignedInterval->onReadStart + align->lastPosition.readPosition;

				rightOfInv->onRefStart = alignedInterval->onRefStart + align->PositionOffset + bestInversionMidpointOnRef;
				rightOfInv->onRefStop = alignedInterval->onRefStart + align->PositionOffset + align->lastPosition.refPosition;
				rightOfInv->isReverse = alignedInterval->isReverse;
			}

		}
	} else {
		verbose(0, true, "Too many peaks detected! %d > %d allowed for read length of %d", checkCount, maxCheckCount, read->length);
		bestResult = SV_NONE;
	}

	return bestResult;
}

int AlignmentBuffer::estimateCorridor(Interval const * interval) {
	int corridor = 8192;
	int onRead = interval->onReadStop - interval->onReadStart;
	int onRef = interval->onRefStop - interval->onRefStart;
	int diff = onRead - onRef;
	int corridorFromDiff = (int) ((abs(diff) * 2.1f));
	int corridorFromLength = (int) ((abs(onRead) * 0.20f));
	corridor = std::min(corridor,
			std::max(corridorFromDiff, corridorFromLength));
	if (pacbioDebug)
		Log.Message("Corridor estimation - onRead: %d, onRef: %d, diff: %d, corridor: %d", onRead, onRef, diff, corridor);

	return corridor;
}

Align * AlignmentBuffer::alignInterval(MappedRead const * const read,
		Interval const * interval, char const * const readSeq,
		size_t const readSeqLen, bool const realign, bool const fullAlignment) {

	Align * align = 0;

	if(readSeq == nullptr) {
		return align;
	}

	Timer alignTimer;
	alignTimer.ST();

	int const minReadSeqLength = 10;
	if (!(llabs(interval->onReadStart - interval->onReadStop) == 0
			|| llabs(interval->onRefStart - interval->onRefStop) == 0
			|| readSeqLen < minReadSeqLength)) {

		int corridor = estimateCorridor(interval);

		int QStart = 0;
		int QEnd = 0;

		if (interval->isReverse) {
			QEnd = interval->onReadStart;
			QStart = read->length - interval->onReadStop;

		} else {
			QStart = interval->onReadStart;
			QEnd = read->length - interval->onReadStop;
		}

		if (pacbioDebug) {
			Log.Message("Computing alignment - Start pos: %d, Length: %d", interval->onReadStart, readSeqLen);
		}
		align = computeAlignment(interval, corridor, readSeq, readSeqLen,
				QStart, QEnd, read->length, read, realign, fullAlignment,
				false);
	} else {
			verbose(0, "Tried to align invalid interval:", interval);
	}
	alignTime += alignTimer.ET();

	return align;
}

unique_ptr<char []> AlignmentBuffer::extractReadSeq(int const readSeqLen,
		int const onReadStart, bool const isReverse, MappedRead* read,
		bool const revComp) {

	// > 500000 very basic check for overflows (this is terrible)
	if(readSeqLen <= 0 || readSeqLen > 200000000) {
		return 0;
	}

	unique_ptr<char[]> readSeq(new char[readSeqLen + 1]);
	if (isReverse) {
		read->computeReverseSeq();
		computeReverseSeq(read->Seq + onReadStart, readSeq.get(), readSeqLen);
		readSeq[readSeqLen] = '\0';
	} else {
		strncpy(readSeq.get(), read->Seq + onReadStart, readSeqLen);
		readSeq[readSeqLen] = '\0';
	}

	if(revComp) {
		unique_ptr<char[]> tmp(new char[readSeqLen + 1]);
		computeReverseSeq(readSeq.get(), tmp.get(), readSeqLen);
		tmp[readSeqLen] = '\0';
		return tmp;
	} else {
		return readSeq;
	}
}


unique_ptr<char []> AlignmentBuffer::extractReadSeq(int const readSeqLen,
		Interval const * interval, MappedRead* read,
		bool const revComp) {
	return extractReadSeq(readSeqLen, interval->onReadStart, interval->isReverse, read, revComp);
}

int AlignmentBuffer::realign(int const svTypeDetected,
		Interval const * interval, Interval * leftOfInv, Interval * rightOfInv,
		MappedRead * read, Align * tmpAling, int & alignIxndex,
		LocationScore * tmpScore, int mq) {

	verbose(0, true, "***********************************");
	verbose(0, true, "Realigning read!");
	verbose(0, "Interval: ", interval);

	float realignScore = 0.0f;
	int svTypeResult = 0;

	/**********************************************/
	//Align part of read that is left of inversion
	/**********************************************/
	verbose(0, "Left Interval:", leftOfInv);
	int readSeqLen = leftOfInv->onReadStop - leftOfInv->onReadStart;

	Align * alignLeft = 0;
	alignLeft = alignInterval(read, leftOfInv, extractReadSeq(readSeqLen, leftOfInv, read).get(), readSeqLen, true, false);

	LocationScore scoreLeft;
	Align * alignRight = 0;
	LocationScore scoreRight;
	Align * alignInv = 0;
	LocationScore scoreInv;



	// Defines the inverted part of the alignment
	// start and end are computed from the alignments
	// left and right to it
	Interval * inv = new Interval();
	if (alignLeft != 0 && alignLeft->Score > 0.0f) {
		alignLeft->MQ = mq;

		scoreLeft.Location.m_Location = leftOfInv->onRefStart
		+ alignLeft->PositionOffset;//- (corridor >> 1); handled in computeAlingment
		scoreLeft.Location.setReverse(leftOfInv->isReverse);
		scoreLeft.Score.f = alignLeft->Score;
		/**********************************************/
		//QEnd from first read part, gives length of
		//inversion that was covered by left part
		/**********************************************/
		// Set start of inverted read part
		inv->onReadStart = read->length - alignLeft->QEnd;
		inv->onRefStart = scoreLeft.Location.m_Location
		+ alignLeft->lastPosition.refPosition;
		inv->isReverse = !leftOfInv->isReverse;

		realignScore += alignLeft->Score;

		/**********************************************/
		//Align part of read that is right of inversion
		/**********************************************/
		verbose(0, "Right Interval: ", rightOfInv);
		readSeqLen = rightOfInv->onReadStop - rightOfInv->onReadStart;
		alignRight = alignInterval(read, rightOfInv, extractReadSeq(readSeqLen, rightOfInv, read).get(), readSeqLen, true, false);

		if (alignRight != 0 && alignRight->Score > 0.0f) {
			alignRight->MQ = mq;

			scoreRight.Location.m_Location = rightOfInv->onRefStart
			+ alignRight->PositionOffset;//- (corridor >> 1); handled in computeAlingment
			scoreRight.Location.setReverse(rightOfInv->isReverse);
			scoreRight.Score.f = alignRight->Score;
			// Set end of inverted read part
			/**********************************************/
			//QStart from second read part, gives length
			//of inversion that was covered by second read
			/**********************************************/
			inv->onReadStop = alignRight->QStart;
			inv->onRefStop = scoreRight.Location.m_Location
			+ alignRight->firstPosition.refPosition;

			realignScore += alignRight->Score;

			/**********************************************/
			//Align inverted part, length was inferred
			//from left and right part
			/**********************************************/
			if(!inv->isReverse) {
				// Reverse coordinates on read
				uloc tmp = read->length - inv->onReadStart;
				inv->onReadStart = read->length - inv->onReadStop;
				inv->onReadStop = tmp;
			}
			// Extend
			int inversionLength = llabs(inv->onRefStop - inv->onRefStart);
			int const minInversionLength = Config.getMinInversionLength();

			verbose(0, true, "Length of inversion is %d", inversionLength);
			verbose(0, "Inverted Interval: ", inv);
			if(inversionLength > minInversionLength) {
				readSeqLen = inv->onReadStop - inv->onReadStart;
				alignInv = alignInterval(read, inv, extractReadSeq(readSeqLen, inv, read).get(), readSeqLen, true, true);

				if(pacbioDebug && alignInv != 0) {
					Log.Message("Inversion alignment score: %f", alignInv->Score);
				}

				Align * alignInvRev = 0;
				alignInvRev = alignInterval(read, inv, extractReadSeq(readSeqLen, inv, read, true).get(), readSeqLen, true, true);

				if(pacbioDebug && alignInvRev != 0) {
					Log.Message("Inversion alignment score: %f", alignInvRev->Score);
				}

				if (alignInv != 0 && alignInv->Score > 0.0f && alignInv->getAlignedReadBp(read->length) > minInversionLength && (alignInvRev == 0 || alignInvRev->Score < alignInv->Score)) {
					alignInv->MQ = mq;
//					alignInv->svType = SV_INVERSION;

					scoreInv.Location.m_Location = inv->onRefStart
					+ alignInv->PositionOffset;	//- (corridor >> 1); handled in computeAlingment
					scoreInv.Location.setReverse(inv->isReverse);
					scoreInv.Score.f = alignInv->Score;

					SequenceLocation inversionStartRef;
					inversionStartRef.m_Location = inv->onRefStart;
					SequenceProvider.convert(inversionStartRef);
					SequenceLocation inversionStopRef;
					inversionStopRef.m_Location = inv->onRefStop;
					SequenceProvider.convert(inversionStopRef);

					realignScore += alignInv->Score;

					svTypeResult = SV_INVERSION;
				} else {
					if (pacbioDebug) {
						Log.Message("Alignment of inverted part failed or score not high enough. Reporting translocation.");
					}
					svTypeResult = SV_TRANSLOCATION;
				}
				if (alignInvRev != 0) {
					alignInvRev->clearBuffer();
					alignInvRev->clearNmPerPosition();
					delete alignInvRev;
					alignInvRev = 0;
				}
			} else {
				if (pacbioDebug) {
					Log.Message("Inversion/Translocation too short.");
				}
				svTypeResult = SV_NONE;
			}
		} else {
			if (pacbioDebug) {
				Log.Message("Alignment right failed");
			}
			svTypeResult = SV_NONE;
		}
	} else {
		if (pacbioDebug) {
			Log.Message("Alignment left failed");
		}

		svTypeResult = SV_NONE;
	}

	delete inv;
	inv = 0;

	if(svTypeResult == SV_NONE || (alignLeft == 0 || alignLeft->Score <= 0.0f)
			|| (alignRight == 0 || alignRight->Score <= 0.0f)) {
		if(alignLeft != 0) {
			alignLeft->clearBuffer();
			alignLeft->clearNmPerPosition();
			delete alignLeft;
			alignLeft = 0;
		}
		if(alignRight != 0) {
			alignRight->clearBuffer();
			alignRight->clearNmPerPosition();
			delete alignRight;
			alignRight = 0;
		}
		if(alignInv != 0) {
			alignInv->clearBuffer();
			alignInv->clearNmPerPosition();
			delete alignInv;
			alignInv = 0;
		}
		svTypeResult = SV_NONE;
	} else {
		// Add
		alignLeft->clearNmPerPosition();
		tmpAling[alignIxndex] = *alignLeft;
		delete alignLeft;
		alignLeft = 0;
		tmpScore[alignIxndex] = scoreLeft;
		tmpAling[alignIxndex].mappedInterval = getIntervalFromAlign(&tmpAling[alignIxndex], &tmpScore[alignIxndex], alignIxndex, read->length);
		alignIxndex += 1;
		read->Calculated += 1;
		alignRight->clearNmPerPosition();
		tmpAling[alignIxndex] = *alignRight;
		delete alignRight;
		alignRight = 0;
		tmpScore[alignIxndex] = scoreRight;
		tmpAling[alignIxndex].mappedInterval = getIntervalFromAlign(&tmpAling[alignIxndex], &tmpScore[alignIxndex], alignIxndex, read->length);
		alignIxndex += 1;
		read->Calculated += 1;

		if(svTypeResult == SV_INVERSION && alignInv != 0) {
			alignInv->clearNmPerPosition();
			tmpAling[alignIxndex] = *alignInv;
			delete alignInv;
			alignInv = 0;
			tmpScore[alignIxndex] = scoreInv;
			tmpAling[alignIxndex].mappedInterval = getIntervalFromAlign(&tmpAling[alignIxndex], &tmpScore[alignIxndex], alignIxndex, read->length);
			alignIxndex += 1;
			read->Calculated += 1;
		} else {
			if(alignInv != 0) {
				alignInv->clearBuffer();
				alignInv->clearNmPerPosition();
				delete alignInv;
				alignInv = 0;
			}
		}
	}

	if(pacbioDebug) {
		Log.Message("Realign score: %f", realignScore);
	}
	return svTypeResult;
}

bool satisfiesConstraints(Align * align, int const readLength) {
	//TODO: check threshold
	static float minResidues = 50.0f;//Config.getMinResidues();
	static float minIdentity = Config.getMinIdentity();

	if (minResidues <= 1.0f) {
		minResidues = readLength * minResidues;
	}
	return (align->Score > 0.0f) && (align->Identity >= minIdentity) && ((float) (readLength - align->QStart - align->QEnd) >= minResidues);
}


void AlignmentBuffer::alignSingleOrMultipleIntervals(MappedRead * read, Interval const * const interval, LocationScore * tmp, Align * tmpAling, int & alignIndex) {

	int readSeqLen = interval->onReadStop - interval->onReadStart;
	auto readPartSeq = extractReadSeq(readSeqLen, interval, read);

	if (readPartSeq != 0) {
		Align * align = alignInterval(read, interval, readPartSeq.get(), readSeqLen, false, false);
		if (align != 0) {
			if (align->Score > 0.0f) {
				int svType = SV_NONE;

				if (Config.getSmallInversionDetection() || Config.getLowQualitySplit()) {
					Interval * leftOfInv = new Interval();
					Interval * rightOfInv = new Interval();
					svType = SV_UNKNOWN;
					bool inversionAligned = false;
					svType = detectMisalignment(align, interval, readPartSeq.get(), leftOfInv, rightOfInv, read);

					if (svType != SV_NONE) {
						int mq = computeMappingQuality(*align, read->length);
						int assumedSvType = svType;
						svType = realign(svType, interval, leftOfInv, rightOfInv, read, tmpAling, alignIndex, tmp, mq);
					} else {
						verbose(0, true, "No SV detected!");
					}
					if (leftOfInv != 0) {
						delete leftOfInv;
						leftOfInv = 0;
					}
					if (rightOfInv != 0) {
						delete rightOfInv;
						rightOfInv = 0;
					}
				}

				if (svType == SV_NONE) {
					verbose(0, true, "No SV was detected. Keeping single alignment!");
					/**********************************************/
					// No inversion detected
					/**********************************************/
					if (satisfiesConstraints(align, read->length)) {
						align->MQ = computeMappingQuality(*align, read->length);
						align->clearNmPerPosition();

						tmpAling[alignIndex] = *align;
						delete align;
						align = 0;

						tmp[alignIndex].Location.m_Location = interval->onRefStart + tmpAling[alignIndex].PositionOffset;					//- (corridor >> 1); handled in computeAlingment
						tmp[alignIndex].Location.setReverse(interval->isReverse);
						tmp[alignIndex].Score.f = tmpAling[alignIndex].Score;

						tmpAling[alignIndex].mappedInterval = getIntervalFromAlign(&tmpAling[alignIndex], &tmp[alignIndex], alignIndex, read->length);

						read->Calculated += 1;
						alignIndex += 1;
					} else {
						align->clearBuffer();
						align->clearNmPerPosition();
						delete align;
						align = 0;
						verbose(0, true, "Alignment did not satisfiesConstraints");
					}
				} //else {
					if (align != 0) {
						align->clearBuffer();
						align->clearNmPerPosition();
						delete align;
						align = 0;
					}
//				}
			} else {
				if (align != 0) {
					align->clearBuffer();
					align->clearNmPerPosition();
					delete align;
					align = 0;
				}

				verbose(0, true, "Alignment failed");
			}
		}
	} else {
		Log.Message("Extracting read sequence failed for read %s. Please report this on https://github.com/philres/ngmlr", read->name);
	}
}

int AlignmentBuffer::computeMappingQuality(Align const & alignment, int readLength) {

	std::vector<IntervalTree::Interval<int> > results;

	verbose(0, true, "Computing mapping quality:");

	readCoordsTree->findOverlapping(alignment.QStart, readLength - alignment.QEnd, results);
	int mqSum = 0;
	int mqCount = 0;
	for (int j = 0; j < results.size(); ++j) {
		//verbose(1, false, "%d, ", results[j].value);
		mqSum += results[j].value;
		mqCount += 1;
	}
//	verbose(1, true, "");
	if (mqCount == 0) return 0;
	verbose(1, true, "%d / %d = %d", mqSum, mqCount, (int) (mqSum * 1.0f / mqCount));
	return (int) (mqSum * 1.0f / mqCount);

//		std::vector<IntervalTree::Interval<int> > results;
//
//	readCoordsTree->findOverlapping(alignment.QStart,
//			readLength - alignment.QEnd, results);
//	int mq = 0;
//
//	verbose(0, true, "Computing mq (readlength %d): ", readLength);
//	int * mqs = new int[results.size()];
//	for (int j = 0; j < results.size(); ++j) {
//		mqs[j] = results[j].value;
//		verbose(1, false, "%d, ", mqs[j]);
//	}
//	std::sort(mqs, mqs + results.size(), std::greater<int>());
//
//	int length = std::min((int)results.size(), std::max(2, (int)(results.size() * 0.2f + 0.5f)));
//	verbose(1, false, "\nUsing (%d): ", length);
//
//	int mqSum = 0;
//	int mqCount = 0;
//	for (int j = 0; j < length; ++j) {
//		mqSum += mqs[j];
//		mqCount += 1;
//		verbose(0, false, "%d, ", mqs[j]);
//	}
//	mq = (int) (mqSum * 1.0f / mqCount);
//	verbose(1, true, "\nMapping quality: %d", mq);
//
//	delete[] mqs; mqs = 0;
//
//	return mq;
}

bool sortMappedSegements(Interval * a,
		Interval * b) {
	return a->score > b->score;
//	int aLen = a.value->onReadStop - a.value->onReadStart;
//	int bLen = b.value->onReadStop - b.value->onReadStart;
//	return (aLen == bLen && a.value->score > b.value->score) || aLen > bLen;
}

bool isContainedOnRead(Interval * shortInterval, Interval * longInterval,
		float minOverlap) {

	int threshold = minOverlap
			* abs(shortInterval->onReadStop - shortInterval->onReadStart);

	int overlap = 0;
	if (longInterval->onReadStart > shortInterval->onReadStart) {
//		Interval * tmp = longInterval;
//		longInterval = shortInterval;
//		shortInterval = tmp;
		overlap = std::max(0,
				shortInterval->onReadStop - longInterval->onReadStart);
	} else {
		overlap = std::max(0,
				longInterval->onReadStop - shortInterval->onReadStart);
	}

//	int overlap = abs(longInterval->onReadStart - shortInterval->onReadStart)
//			+ abs(longInterval->onReadStop - shortInterval->onReadStop);
	return overlap > threshold;
}

bool isFullyContainedOnRead(Interval * shortInterval, Interval * longInterval) {

//	int threshold = 0.01f * (longInterval->onReadStop - longInterval->onReadStart);
	int threshold = 10;
	return (longInterval->onReadStart - shortInterval->onReadStart) <= threshold
			&& (longInterval->onReadStop - shortInterval->onReadStop)
					>= -threshold;
}

bool isFullyContainedOnRef(Interval * shortInterval, Interval * longInterval) {

//	int threshold = 0.10f * (longInterval->onReadStop - longInterval->onReadStart);
	int threshold = 10;
	return (longInterval->onRefStart - shortInterval->onRefStart) <= threshold
			&& (longInterval->onRefStop - shortInterval->onRefStop)
					>= -threshold;
}

bool isFullyContained(Interval * shortInterval, Interval * longInterval) {

	return isFullyContainedOnRead(shortInterval, longInterval)
			&& isFullyContainedOnRef(shortInterval, longInterval);
}

bool isValidOverlap(Interval * shortInterval, Interval * longInterval) {
	return true;
}

bool isValidOverlapRef(Interval * shortInterval, Interval * longInterval) {
//	return shortInterval->isReverse == longInterval->isReverse
//			&& shortInterval->onRefStart <= longInterval->onRefStop
//			&& shortInterval->onRefStop >= longInterval->onRefStart;
	return shortInterval->isReverse == longInterval->isReverse
			&& (longInterval->onRefStop - shortInterval->onRefStart) > 50
			&& (shortInterval->onRefStop - longInterval->onRefStart) > 50;
}

bool isValidOverlapRead(Interval * shortInterval, Interval * longInterval) {
//	return shortInterval->isReverse == longInterval->isReverse
//			&& shortInterval->onRefStart <= longInterval->onRefStop
//			&& shortInterval->onRefStop >= longInterval->onRefStart;
	return shortInterval->isReverse == longInterval->isReverse
			&& (longInterval->onReadStop - shortInterval->onReadStart) > 50
			&& (shortInterval->onReadStop - longInterval->onReadStart) > 50;
}

float getBestSegmentCombination(int const maxLength, std::vector<Interval *> const & mappedSegements, std::vector<int> & segmentList) {

	float result = 0.0f;

	// field i = score for best fragment (combination) in range 0 to i
	float * bestScore = new float[maxLength];
	int * lastBest = new int[maxLength];
	int * lastFragment = new int[maxLength];

	int const maxOverlap = 50;

	bestScore[0] = 0;
	lastFragment[0] = -1;
	lastBest[0] = 0;
	for (int i = 1; i < maxLength; ++i) {

		// Init
		bestScore[i] = bestScore[i - 1];
		lastFragment[i] = lastFragment[i - 1];
		lastBest[i] = lastBest[i - 1];

		// Find fragment (combination) with the highest score that fits into interval 0 to i
		for (int j = 0; j < mappedSegements.size(); ++j) {

			if (!mappedSegements[j]->isProcessed && mappedSegements[j]->onReadStop <= i && abs(mappedSegements[j]->onReadStop - mappedSegements[j]->onReadStart) > maxOverlap) {

				int const start = std::min(maxLength, mappedSegements[j]->onReadStart + maxOverlap);
				float currentScore = mappedSegements[j]->score + bestScore[start];
				if (currentScore > bestScore[i]) {
					bestScore[i] = currentScore;
					lastFragment[i] = j;
					lastBest[i] = start;
				}
			} else {
				// Fragment doesn't fit into range 0 to i
			}

		}

	}

	int i = maxLength - 1;
	result = bestScore[i];

	while (lastFragment[i] > -1) {

		//Log.Message("Fragment: %d", lastFragment[i]);
		segmentList.push_back(lastFragment[i]);
		i = lastBest[i];
	}

	delete[] bestScore;
	bestScore = 0;
	delete[] lastBest;
	lastBest = 0;
	delete[] lastFragment;
	lastFragment = 0;

	return result;
}

Interval * AlignmentBuffer::getIntervalFromAlign(Align const * const align, LocationScore const * const score, int const i, int const readLength) {
	int diffOnRef = align->lastPosition.refPosition - align->firstPosition.refPosition;

	Interval * mappedSegement = new Interval();
	mappedSegement->id = i;
	mappedSegement->onRefStart = score->Location.m_Location;
	mappedSegement->onRefStop = score->Location.m_Location + diffOnRef;
	mappedSegement->isReverse = score->Location.isReverse();
	mappedSegement->isProcessed = false;
	mappedSegement->score = align->Score;
	if (mappedSegement->isReverse) {
		// If mapped to the reverse strand, QStart and QEnd have to be translated to
		// plus strand
		mappedSegement->onReadStart = align->QEnd;
		mappedSegement->onReadStop = readLength - align->QStart - 1;
	} else {
		mappedSegement->onReadStart = align->QStart;
		mappedSegement->onReadStop = readLength - align->QEnd - 1;
	}

	return mappedSegement;
}

bool AlignmentBuffer::reconcileRead(ReadGroup * group) {

	bool mapped = false;

	std::vector<Interval *> mappedSegements;

	MappedRead * read = group->fullRead;

	// Not the best strategy, however only used for debug
	// plots at the moment (ngm itself doesn't rely on read strand information)
	bool readIsReverse = read->Scores[0].Location.isReverse();

	if (pacbioDebug) {
		Log.Message("Reconciling read: %s with length %d (%d)", read->name, read->length, read->ReadId);
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
		mappedSegements.push_back(mappedSegement);

		if (stdoutPrintMappedSegments) {
			fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", read->name, i,
					read->length, mappedSegement->onReadStart,
					mappedSegement->onReadStop,
					read->Scores[i].Location.isReverse(),
					read->Alignments[i].MQ, read->Alignments[i].Score);
		}
	}


//	// Sort segments by score
//	std::sort(mappedSegements.begin(), mappedSegements.end(),
//			sortMappedSegements);


	// DEBUG code: print segments for dotplot and log
	for (int i = 0; i < mappedSegements.size(); ++i) {
		if (pacbioDebug) {
			Log.Message("%d: Read: %d - %d, Ref: %llu - %llu (%d) - %f (%f, %d) - %d", mappedSegements[i]->id, mappedSegements[i]->onReadStart, mappedSegements[i]->onReadStop, mappedSegements[i]->onRefStart, mappedSegements[i]->onRefStop, mappedSegements[i]->isReverse, mappedSegements[i]->score, read->Alignments[mappedSegements[i]->id].Identity, (read->length - read->Alignments[mappedSegements[i]->id].QStart - read->Alignments[mappedSegements[i]->id].QEnd), mappedSegements[i]->isProcessed);
		}

		if(!readIsReverse) {
			printDotPlotLine(read->ReadId, read->name,
					mappedSegements[i]->onReadStart,
					mappedSegements[i]->onReadStop,
					mappedSegements[i]->onRefStart,
					mappedSegements[i]->onRefStop,
					mappedSegements[i]->score,
					mappedSegements[i]->isReverse,
					DP_TYPE_RESULT + mappedSegements[i]->id, DP_STATUS_OK);
		} else {
			printDotPlotLine(read->ReadId, read->name,
					mappedSegements[i]->onReadStop,
					mappedSegements[i]->onReadStart,
					mappedSegements[i]->onRefStart,
					mappedSegements[i]->onRefStop,
					mappedSegements[i]->score,
					mappedSegements[i]->isReverse,
					DP_TYPE_RESULT + mappedSegements[i]->id, DP_STATUS_OK);
		}
	}


	int const maxLength = read->length;


	std::vector<int> bestSegments;
	float bestScore = getBestSegmentCombination(maxLength, mappedSegements, bestSegments);

	if (pacbioDebug) {
		Log.Message("Best score: %f", bestScore);
	}

	float topFragmentScore = 0.0f;
	int fragmentWithHighestScore = 0;

	int alignedBpSum = 0;
	for(int i = 0; i < bestSegments.size(); ++i) {
		int index = bestSegments[i];

		if (pacbioDebug) {
			Log.Message("\tFragment: %d", index);
		}
		mappedSegements[index]->isProcessed = true;

		alignedBpSum += (mappedSegements[index]->onReadStop - mappedSegements[index]->onReadStart);
		if(mappedSegements[index]->score > topFragmentScore) {
			fragmentWithHighestScore = index;
			topFragmentScore = mappedSegements[index]->score;
		}
	}
	// Set the alignment with the highest score as primary
	if(bestSegments.size() > 0) {
		read->Alignments[mappedSegements[fragmentWithHighestScore]->id].primary = true;
	}
	float aligned = alignedBpSum * 1.0f / read->length;
	if(pacbioDebug) {
		Log.Message("Aligned %.2f%% of read", aligned * 100.0f);
	}
	mapped = (Config.getMinResidues() < 1.0) ? aligned > Config.getMinResidues() : alignedBpSum > Config.getMinResidues();

	if (pacbioDebug) {
		Log.Message("%f > %f = %d", aligned, Config.getMinResidues(), mapped);
	}
	NGM.Stats->avgAlignPerc += aligned;
	bestSegments.clear();
	float secondBestScore = getBestSegmentCombination(maxLength, mappedSegements, bestSegments);
	if(pacbioDebug) {
		Log.Message("Second best score: %f", secondBestScore);
	}

	/**
	 * Filter short isolated intervals
	 * TODO: improve
	 */
	verbose(0, true, "Filter short isolated alignments");
	int const minOnReadLength = 1000;

	for (int i = 0; i < mappedSegements.size(); ++i) {
		Interval * a = mappedSegements[i];
		if (a->isProcessed) {
			verbose(1, "Testee: ", a);
			bool keep = a->lengthOnRead() > std::min(minOnReadLength, (int)(read->length * 0.5f));
			for (int j = 0; j < mappedSegements.size() && !keep; ++j) {
				Interval * b = mappedSegements[j];
				if (b != 0 && b->isProcessed) {
					verbose(2, "Compare: ", b);
					int distance = getDistanceOnRead(a, b);
					loc distanceOnRef = (b->onRefStart < a->onRefStart) ? std::max(0ll, a->onRefStart - b->onRefStop) : std::max(0ll, b->onRefStart - a->onRefStop);
					loc const maxDistance = a->lengthOnRead();
					keep = (distance < maxDistance || distanceOnRef < maxDistance) && b->lengthOnRead() > std::min(minOnReadLength, (int)(read->length * 0.5f));
					verbose(2, true, "dist: %d, a length * 2: %d, b length: %d -> %d", distance, a->lengthOnRead() * 2, b->lengthOnRead(), keep);
				}
			}
			if (!keep) {
				verbose(2, true, "Delete");
				a->isProcessed = false;
			}
		}
	}

	// Mark remaining as skipped (will not be written to SAM file)
	for (int i = 0; i < mappedSegements.size(); ++i) {
		if(!mappedSegements[i]->isProcessed) {
			read->Alignments[mappedSegements[i]->id].skip = true;
		}
	}

	int segmentCount = 0;
	// Set 0x1 (full read is aligned) + debug code
	for (int i = 0; i < read->Calculated; ++i) {
		if (!read->Alignments[mappedSegements[i]->id].skip) {

			if(aligned > 0.95f) {
				read->Alignments[mappedSegements[i]->id].setBitFlag(0x2);
			}

			segmentCount += 1;
			if (!readIsReverse) {
				printDotPlotLine(read->ReadId, read->name,
						mappedSegements[i]->onReadStart,
						mappedSegements[i]->onReadStop,
						mappedSegements[i]->onRefStart,
						mappedSegements[i]->onRefStop,
						mappedSegements[i]->score,
						mappedSegements[i]->isReverse,
						DP_TYPE_RESULT_CONS + mappedSegements[i]->id,
						DP_STATUS_OK);
			} else {
				printDotPlotLine(read->ReadId, read->name,
						mappedSegements[i]->onReadStop,
						mappedSegements[i]->onReadStart,
						mappedSegements[i]->onRefStart,
						mappedSegements[i]->onRefStop,
						mappedSegements[i]->score,
						mappedSegements[i]->isReverse,
						DP_TYPE_RESULT_CONS + mappedSegements[i]->id,
						DP_STATUS_OK);
			}
		}
	}

	int maxSplits = Config.getMaxSegmentNumberPerKb(read->length);
	mapped = mapped && (segmentCount - 1) <= maxSplits;
	verbose(0, true, "Splits detected: %d < %d = %d (%d)", segmentCount - 1, maxSplits, mapped, read->length);

	// Delete all mapped segments
	for (int i = 0; i < mappedSegements.size(); ++i) {
		delete mappedSegements[i];
		mappedSegements[i] = 0;
	}

	if(mapped) {
		verbose(0, true, "%s (%d) mapped %f with %d fragments", read->name, read->length, aligned * 100.0f, segmentCount);
	}
	return mapped;
}

void sortRead(ReadGroup * group) {

	MappedRead * read = group->fullRead;

	float highestScore = 0.0f;
	int highestScoreIndex = 0;
	for (int i = 0; i < read->Calculated; ++i) {
		if (read->Alignments[i].Score > highestScore) {
			highestScore = read->Alignments[i].Score;
			highestScoreIndex = i;
		}
	}
	if (highestScoreIndex != 0) {
		Align tmp = read->Alignments[0];
		read->Alignments[0] = read->Alignments[highestScoreIndex];
		read->Alignments[highestScoreIndex] = tmp;

		LocationScore tmpScore = read->Scores[0];
		read->Scores[0] = read->Scores[highestScoreIndex];
		read->Scores[highestScoreIndex] = tmpScore;
	}
}

int AlignmentBuffer::getOverlapOnRead(Interval const * a, Interval const * b) {
	return std::max(0, std::min(a->onReadStop, b->onReadStop) - std::max(a->onReadStart, b->onReadStart));
}

int AlignmentBuffer::getDistanceOnRead(Interval const * a, Interval const * b) {
	if(b->onReadStart < a->onReadStart) {
		return std::max(0, a->onReadStart - b->onReadStop);
	} else {
		return std::max(0, b->onReadStart - a->onReadStop);
	}
}

loc AlignmentBuffer::getDistanceOnRef(Interval const * a, Interval const * b) {
	if (b->isReverse) {
		if (b->onRefStop < a->onRefStop) {
			return std::max(0ll, a->onRefStop - b->onRefStart);
		} else {
			return std::max(0ll, b->onRefStop - a->onRefStart);
		}
	} else {
		if (b->onRefStart < a->onRefStart) {
			return std::max(0ll, a->onRefStart - b->onRefStop);
		} else {
			return std::max(0ll, b->onRefStart - a->onRefStop);
		}
	}
}


void AlignmentBuffer::verbose(int const tabs, char const * const s, Interval const * const interval) {
	verbose(tabs, false, "%s", s);
	verbose(0, true, "Interval: %d - %d on read, %lu - %lu on ref, Reverse: %d, Score: %f, Anchors: %d", interval->onReadStart, interval->onReadStop, interval->onRefStart, interval->onRefStop, interval->isReverse, interval->score, interval->anchorLength);
}


void AlignmentBuffer::verbose(int const tabs, bool const newLine, char const * const s, ...) {
	if(pacbioDebug) {
		for(int i = 0; i < tabs; ++i) {
			fprintf(stderr, "\t");
		}
		va_list args;

		va_start(args, s);
		vfprintf(stderr, s, args);
		va_end(args);

		if(newLine) {
			fprintf(stderr, "\n");
		}
	}
}

bool AlignmentBuffer::extendIntervalStop(Interval * interval, int const readBp, int const readLength) {
	bool extended = false;

	verbose(0, true, "## Extending interval stop:");
	verbose(0, "## ", interval);

	_SequenceProvider::Chromosome chr = SequenceProvider.getChrBorders(interval->onRefStart, interval->onRefStop);


	//verbose(0, true, "extendIntervalEnd fro %llu - %llu by %d -- Located on chr %llu - %llu", interval->onRefStart, interval->onRefStop, readBp, chr.start, chr.end);

	if (chr.start != 0 || chr.end != 0) {
		verbose(0, true, "extendIntervalStop - Located on chr %llu %llu", chr.start, chr.end);

		double lengthRatio = std::min(1.0f, interval->lengthOnRead() * 1.0f / interval->lengthOnRef() * 1.0f);

		int extendOnRead = std::min(readLength - interval->onReadStop, readBp);
		int extendOnRef = (int) round(extendOnRead / lengthRatio);

		loc maxExtendOnRef = interval->onRefStop > chr.end ? 0 : chr.end - interval->onRefStop;
		if (interval->isReverse) {
			maxExtendOnRef = interval->onRefStop < chr.start ? 0 : interval->onRefStop - chr.start;
		}

		if (extendOnRef > maxExtendOnRef) {
			extendOnRef = maxExtendOnRef;
			extendOnRead = std::min(extendOnRead, std::max(0, (int) round(extendOnRef * lengthRatio) - 1));
		}

		verbose(1, true, "Min/Max extend on ref: %d/%lld", extendOnRef, maxExtendOnRef);
		interval->onReadStop += extendOnRead;
		extended = true;
		if (interval->isReverse) {
			interval->onRefStop -= extendOnRef;
		} else {
			interval->onRefStop += extendOnRef;
		}
	} else {
		verbose(0, true, "Could not get chromosome for interval");
		getchar();
	}

	return extended;
}

bool AlignmentBuffer::extendIntervalStart(Interval * interval, int const readBp) {
	bool extended = false;

	_SequenceProvider::Chromosome chr = SequenceProvider.getChrBorders(interval->onRefStart, interval->onRefStop);

	verbose(0, true, "extendIntervalStart fro %llu - %llu by %d -- Located on chr %llu - %llu", interval->onRefStart, interval->onRefStop, readBp, chr.start, chr.end);

	if (chr.start != 0 || chr.end != 0) {

		double lengthRatio = std::min(1.0f, interval->lengthOnRead() * 1.0f / interval->lengthOnRef() * 1.0f);

		int extendOnRead = std::min(interval->onReadStart, readBp);
		int extendOnRef = (int) round(extendOnRead / lengthRatio);

		loc maxExtendOnRef = interval->onRefStart < chr.start ? 0 : interval->onRefStart - chr.start;
		if (interval->isReverse) {
			maxExtendOnRef = interval->onRefStart > chr.end ? 0 : chr.end - interval->onRefStart;
		}
		verbose(1, true, "Min/Max extend on ref: %d/%lld", extendOnRef, maxExtendOnRef);
		if (extendOnRef > maxExtendOnRef) {
			extendOnRef = maxExtendOnRef;
			extendOnRead = std::min(extendOnRead, std::max(0, (int) round(extendOnRef * lengthRatio) - 1));
		}

		interval->onReadStart -= extendOnRead;
		extended = true;
		if (interval->isReverse) {
			interval->onRefStart += extendOnRef;
		} else {
			interval->onRefStart -= extendOnRef;
		}
	} else {
		verbose(0, true, "Could not get chromosome for interval");
//		getchar();
	}

	return extended;
}


bool AlignmentBuffer::shortenIntervalStart(Interval * interval, int const readBp) {

	bool shortened = false;

	if (interval->onReadStart < interval->onReadStop) {
		double lengthRatio = std::max(1.1f, interval->lengthOnRead() * 1.0f / interval->lengthOnRef() * 1.0f);
//		lengthRatio = 1.1;

		int refBp = (int) round(readBp / lengthRatio);

		if (readBp < interval->lengthOnRead() && refBp < interval->lengthOnRef()) {

			interval->onReadStart = interval->onReadStart + readBp;
			interval->onRefStart = interval->isReverse ? interval->onRefStart - refBp : interval->onRefStart + refBp;
			shortened = true;
		}
	}

	return shortened;
}

bool AlignmentBuffer::shortenIntervalEnd(Interval * interval, int const readBp) {

	bool shortened = false;

	if (interval->onReadStart < interval->onReadStop) {
		double lengthRatio = std::max(1.1f, interval->lengthOnRead() * 1.0f / interval->lengthOnRef() * 1.0f);
//		lengthRatio = 1.1;

		int refBp = (int) round(readBp / lengthRatio);

		verbose(1, true, "%llu %d", interval->onRefStop, refBp);

		if (readBp < interval->lengthOnRead() && refBp < interval->lengthOnRef()) {

			interval->onReadStop = interval->onReadStop - readBp;
			interval->onRefStop = interval->isReverse ? interval->onRefStop + refBp : interval->onRefStop - refBp;
			shortened = true;
		}
	}

	return shortened;
}

float AlignmentBuffer::scoreInterval(Interval * interval, MappedRead * read) {
	float score = 0.0f;

	int const tabs = 3;
	if (interval->onReadStart < interval->onReadStop) {
		int onReadStart = interval->onReadStart;
		int onReadStop = interval->onReadStop;

		auto readSeq = extractReadSeq(interval->lengthOnRead(), onReadStart, interval->isReverse, read, false);
		if (readSeq != nullptr) {
			verbose(0, "", interval);
			loc onRefStart = interval->isReverse ? interval->onRefStop : interval->onRefStart;
			loc onRefStop = interval->isReverse ? interval->onRefStart : interval->onRefStop;

			if (onRefStart < onRefStop) {
				int refSeqLength = 0;
				auto refSeq = extractReferenceSequenceForAlignment(onRefStart, onRefStop, refSeqLength);

				if (readSeq != nullptr) {
					if (refSeq != nullptr) {
						verbose(tabs, true, "RefSeq: %s", refSeq);
						verbose(tabs, true, "ReadSeq: %s", readSeq.get());

						overlapCheckAligner->SingleScore(0, 0, refSeq, readSeq.get(), score, 0);

						delete[] refSeq;
						refSeq = 0;
					}
				}
			}
		}
	}
	return score;
}

void AlignmentBuffer::processShortRead(MappedRead * read) {
	// Read is shorter then anchor length

	verbose(0, true, "Read too short for long read cLIS");
	verbose(0, true, "Processing ShortRead: %d - %s (lenght %d)", read->ReadId, read->name, read->length);

	if (read->numScores() > 0) {

		//part->Scores[0].Score.f;

		read->Alignments = new Align[read->numScores()];
		LocationScore * tmpScore = new LocationScore[read->numScores()];
		int alignIndex = 0;

		int lastScore = 0;
		for (int k = 0; k < read->numScores() && ((int) read->Scores[k].Score.f >= lastScore || alignIndex < 2); ++k) {
			LocationScore loc = read->Scores[k];

			lastScore = read->Scores[k].Score.f;

			if (pacbioDebug) {
				Log.Message("Hit found with score %f", loc.Score.f);
			}

			Interval * interval = new Interval();

			int const refExtend = read->length * 0.15f;

			interval->onReadStart = 0;
			interval->onReadStop = read->length;
			interval->onRefStart = loc.Location.m_Location
			- refExtend;
			interval->onRefStop = loc.Location.m_Location + read->length + refExtend;
			interval->isReverse = loc.Location.isReverse();

			verbose(0, "Aligning interval: ", interval);

			int const shortReadCorridor = readPartLength * 1 + 2 * refExtend;

			Align * align = 0;
			auto readPartSeq = extractReadSeq(read->length, interval, read);

//			Align * align = alignInterval(read, interval, read->Seq, read->length, false);
			if(readPartSeq != nullptr) {
				align = computeAlignment(interval, shortReadCorridor, readPartSeq.get(),
						read->length, 0, 0, read->length, read, false, false, true);

			}


			bool mapped = align != 0 && align->Score > 0.0f;
			if(mapped) {
				if (Config.getMinResidues() < 1.0) {
					mapped = ((read->length - align->QStart - align->QEnd) * 1.0f / read->length) >
							 Config.getMinResidues();
				} else {
					mapped = (read->length - align->QStart - align->QEnd) > Config.getMinResidues();
				}
			}
			if (mapped) {
				align->clearNmPerPosition();

				align->MQ = read->mappingQlty;
				loc.Location.m_Location = interval->onRefStart
				+ align->PositionOffset;//- (corridor >> 1); handled in computeAlingment
				loc.Location.setReverse(interval->isReverse);
				loc.Score.f = align->Score;

				read->Alignments[alignIndex] = *align;

				read->Calculated += 1;
				delete align;
				align = 0;
				tmpScore[alignIndex] = loc;
				alignIndex += 1;
			} else {
				verbose(0, true, "Alignment failed");
				if(align != 0) {
					align->clearBuffer();
					align->clearNmPerPosition();
					delete align;
					align = 0;
				}

			}

			delete interval;
			interval = 0;

		}

		if (read->Scores != 0) {
			delete[] read->Scores;
			read->Scores = 0;
		}
		read->AllocScores(tmpScore, alignIndex);
		read->Calculated = alignIndex;

		delete[] tmpScore;
		tmpScore = 0;

		if (alignIndex > 0) {
			read->Alignments[0].primary = true;
			WriteRead(read, true);
		} else {
			WriteRead(read, false);
		}
	} else {
		WriteRead(read, false);
	}
}

bool AlignmentBuffer::gapOverlapsWithInterval(Interval * first, Interval * second, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read) {


	Interval gap;

	gap.onReadStart = first->onReadStop + 1;
	gap.onReadStop = std::max(0, second->onReadStart - 1);

	gap.onRefStart = first->onRefStop;
	//gap.onRefStart = first->isReverse ? first->onRefStart : first->onRefStop;
	gap.onRefStop = second->onRefStart;
	//gap.onRefStop = first->isReverse ? second->onRefStop : second->onRefStart;

	gap.isReverse = first->isReverse;

//	verbose(0, "First: ", first);
//	verbose(0, "Second: ", second);
//	verbose(0, "Gap: ", &gap);

	return gapOverlapsWithInterval(&gap, intervalsTree, read);

}

bool AlignmentBuffer::gapOverlapsWithInterval(Interval * gap, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read) {
	float const minOverlap = 50.0;
	int const maxLengthAlignmentCheck = 1000;
	int const minGapLength = (int) (readPartLength * 1.5f);
	bool overlaps = false;

	if (gap->onReadStart < gap->onReadStop) {
		verbose(2, true, "Looking for overlap on read: %d - %d", gap->onReadStart, gap->onReadStop);

		if (gap->lengthOnRead() > minGapLength) {
			std::vector<IntervalTree::Interval<Interval *> > results;
			intervalsTree->findOverlapping(gap->onReadStart, gap->onReadStop, results);

			for (auto node : results) {
				if (!node.value->isProcessed) {
					verbose(3, "", node.value);
					Interval * interval = new Interval();

					interval->onReadStart = gap->onReadStart;
					interval->onReadStop = gap->onReadStop;
					interval->onRefStart = node.value->onRefStart;
					interval->onRefStop = node.value->onRefStop;
					interval->isReverse = node.value->isReverse;

					/**
					 * Long matches to different region (e.g. secondary mappings) should not count as overlap
					 */
					//if (node.value->lengthOnRead() < (2 * readPartLength + gap->lengthOnRead())) {
					if (node.value->lengthOnRead() < ((int)(4.5f * readPartLength) + gap->lengthOnRead())) {
						int overlap = getOverlapOnRead(node.value, gap);
						float overlapPercent = overlap * 100.0f / gap->lengthOnRead();
						verbose(3, "DEBUG: ", interval);
						verbose(3, true, "Overlap: %d (%f %%)", overlap, overlapPercent);

						bool betterScore = true;
						if (overlapPercent > minOverlap) {
							if (read != 0 && gap->lengthOnRead() < maxLengthAlignmentCheck) {
								// Check scores
								float score1 = scoreInterval(interval, read) / interval->lengthOnRead();
								float score2 = scoreInterval(gap, read) / gap->lengthOnRead();
								verbose(3, "true", "Short overlap found. Verifying with alignment score per position: %f vs %f", score1, score2);
								betterScore = score1 > score2;
							}
						}
						overlaps = overlaps || (overlapPercent > minOverlap && betterScore);
					} else {
						verbose(3, true, "Overlapping interval too long (%d >= %d)", node.value->lengthOnRead(), (2 * readPartLength + gap->lengthOnRead()));
					}
					delete interval;
					interval = 0;
				}
			}
		} else {
			verbose(2, true, "Gap too short (%d bp), skipping check", gap->lengthOnRead());
		}
	}
	return overlaps;
}

bool AlignmentBuffer::gapToEndOverlapsWithInterval(Interval * second, int const readLength, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read) {
	Interval gap;

	gap.onReadStart = std::min(readLength, second->onReadStop + 1);
	gap.onReadStop = readLength;

	//	gap.onRefStart = second->isReverse ? first->onRefStart : first->onRefStop;
	//	gap.onRefStop = second->isReverse ? second->onRefStop : second->onRefStart;

	return gapOverlapsWithInterval(&gap, intervalsTree, 0);

}

bool AlignmentBuffer::gapFromStartOverlapsWithInterval(Interval * second, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read) {

	Interval gap;

	gap.onReadStart = 0;
	gap.onReadStop = std::max(0, second->onReadStart - 1);

//	gap.onRefStart = second->isReverse ? first->onRefStart : first->onRefStop;
//	gap.onRefStop = second->isReverse ? second->onRefStop : second->onRefStart;

	return gapOverlapsWithInterval(&gap, intervalsTree, 0);
}

void AlignmentBuffer::closeGapOnRead(Interval * first, Interval * second, int const readLength) {
	if (first->onReadStop < second->onReadStop) {

		int distance = getDistanceOnRead(first, second);

		int maxDistance = (int)(0.25f * readLength); //2 * std::max(first->lengthOnRead(), second->lengthOnRead());

		if (distance > 0 && distance < maxDistance) {
			verbose(3, true, "Closing gap of %d bp", distance);
			verbose(3, "First:  ", first);
			verbose(3, "Second: ", second);

			extendIntervalStop(first, distance, readLength);
			extendIntervalStart(second, distance);

			verbose(3, "New first:  ", first);
			verbose(3, "New second: ", second);
		} else {
			verbose(3, true, "Not closing gap distance > %d", maxDistance);
		}
	}
}

void AlignmentBuffer::extendToReadStart(Interval * interval, int const readLength, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read) {
	//TODO: add interval tree check
	float const maxReadStartExtend = 0.25f;
	int const maxExtend = std::min((int) round(readLength * maxReadStartExtend), interval->lengthOnRead());
	int extend = interval->onReadStart;

	if (extend > 0) {
		if (extend > (readPartLength)) {
			verbose(2, true, "Extending to read start: %d < %d (%d, %d)", extend, maxExtend, (int) round(readLength * maxReadStartExtend), interval->lengthOnRead());
			if (extend <= maxExtend) {
				if (!gapFromStartOverlapsWithInterval(interval, intervalsTree, read)) {
					verbose(2, true, "Extending first interval to read start");
					extendIntervalStart(interval, extend);
				} else {
					verbose(2, true, "Not extending first interval to read start. Overlaps with other interval.");
				}
			} else {
				verbose(2, true, "Not extending first interval to read start. Distance to big.");
			}
		} else {
			verbose(2, true, "Extending first interval to read start");
			extendIntervalStart(interval, extend);
		}
	}
}

void AlignmentBuffer::extendToReadStop(Interval * interval, int const readLength, IntervalTree::IntervalTree<Interval *> * intervalsTree, MappedRead * read) {
	//TODO: add interval tree check
	float const maxReadStopExtend = 0.25f;
	int const maxExtend = std::min((int) round(readLength * maxReadStopExtend), interval->lengthOnRead());
	int extend = readLength - interval->onReadStop;

	if (extend > 0) {
		if (extend > (readPartLength)) {
			verbose(2, true, "Extending to read stop: %d < %d (%d, %d)", extend, maxExtend, (int) round(readLength * maxReadStopExtend), interval->lengthOnRead());
			if (extend <= maxExtend) {
				if (!gapToEndOverlapsWithInterval(interval, readLength, intervalsTree, read)) {
					verbose(2, true, "Extending last interval to read stop");
					extendIntervalStop(interval, extend, readLength);
				} else {
					verbose(2, true, "Not extending last interval to read stop. Overlaps with other interval.");
				}
			} else {
				verbose(2, true, "Not extending last interval to read stop. Distance to big.");
			}
		} else {
			verbose(2, true, "Extending last interval to read stop");
			extendIntervalStart(interval, extend);
		}
	}
}

void AlignmentBuffer::processLongReadLIS(ReadGroup * group) {

	MappedRead * read = group->fullRead;

	int maxAnchorNumber = 10000;
	int const maxNumScores = 1000;

	/**
	 * Get all mapping positions from anchors (non overlapping 256bp parts of reads)
	 * and convert them to Anchor objects
	 * TODO: remove fixed length of 100
	 */
	int maxCandidatesPerReadPart = 100;

	Timer overallTmr;
	overallTmr.ST();

	if (read->length <= readPartLength) {
		Log.Message("Should not be here. Short read in long read pipeline.");
		WriteRead(group->fullRead, false);
	}

	/*
	 * Parts of read that were aligned plus MQ of subalignments
	 */
	std::vector<IntervalTree::Interval<int> > treeIntervals;

//	Log.Message("Processing LongReadLIS: %d - %s (length %d)", read->ReadId, read->name, read->length);

	verbose(0, true, "Processing LongReadLIS: %d - %s (length %d)", read->ReadId, read->name, read->length);
	verbose(0, true, "Anchor list:");

	Anchor * anchorsFwd = new Anchor[maxAnchorNumber];
	int anchorFwdIndex = 0;


//	float * bestScores = new float[group->readNumber];
//	int bestScoresIndex = 0;
//	for (int j = 0; j < group->readNumber; ++j) {
//		MappedRead * part = group->reads[j];
//		if(part->numScores() > 0) {
//			bestScores[bestScoresIndex++] = part->Scores[0].Score.f;
//		}
//	}
//
//	std::sort(bestScores, bestScores + bestScoresIndex);
//
//	for(int i = 0; i < bestScoresIndex; ++i) {
//		verbose(1, false, "%f, ", bestScores[i]);
//	}
//	verbose(0, true, "");
//	float const minScore = bestScores[(int)(bestScoresIndex * 0.8f)];
//	verbose(1, true, "Min score: %f", minScore);
//
//	delete[] bestScores;
//	bestScores = 0;

//	Log.Message("Read: %s", group->fullRead->name);
//	Log.Message("Parts: %d", group->readNumber);
//	for (int j = 0; j < group->readNumber; ++j) {
//		MappedRead * part = group->reads[j];
//
//		int positionOnRead = j * readPartLength;
//
//		Log.Message("\tPart %d:", j);
//		Log.Message("\t\tName: %s", part->name);
//		Log.Message("\t\tMappings: %d", part->numScores());
//		Log.Message("\t\tQuality: %d", part->mappingQlty);
//
//		for (int k = 0; k < part->numScores(); ++k) {
//			Log.Message("\t\t\t%d, %llu, %f", positionOnRead, part->Scores[k].Location.m_Location, part->Scores[k].Score.f );
//		}
//
//	}


	/**
	 * Collect sub-read mapping locations and
	 * build interval tree with read locations
	 * plus mapping quality
	 */
	for (int j = 0; j < group->readNumber; ++j) {
		MappedRead * part = group->reads[j];

		int positionOnRead = j * readPartLength;

		verbose(1, false, "%d: ", positionOnRead);

		if (part->numScores() < maxNumScores) {
			if (part->numScores() > 0) {

				/**
				 * Add mapped read parts + mapping quality.
				 * Used to estimate MQ for read(part)s
				 */
				treeIntervals.push_back(IntervalTree::Interval<int>(positionOnRead, positionOnRead + readPartLength, part->mappingQlty));

				for (int k = 0; k < part->numScores(); ++k) {

					if (stdoutPrintScores) {
						printf("%f\n", part->Scores[k].Score.f);
					}

					/**
					 * Anchor is valid and will be used
					 */
					if (anchorFwdIndex >= (maxAnchorNumber - 1)) {
						verbose(0, true, "Anchor array too small - reallocating.");
						maxAnchorNumber = maxAnchorNumber * 2;
						Anchor * anchorsTmp = new Anchor[maxAnchorNumber];
						memcpy(anchorsTmp, anchorsFwd, anchorFwdIndex * sizeof(Anchor));
						delete[] anchorsFwd;
						anchorsFwd = anchorsTmp;
					}
					Anchor & anchor = anchorsFwd[anchorFwdIndex++];

					anchor.score = part->Scores[k].Score.f;
					anchor.isReverse = part->Scores[k].Location.isReverse();
					anchor.type = DP_STATUS_OK;
					anchor.isUnique = part->numScores() == 1; // || anchor.score > minScore;

					/**
					 * It would be best to convert reads or read parts that map to the negative strand
					 * Immediately to plus strand.
					 * If somebody tells me how to do this properly, I'll happily change this.
					 * For now Anchors can be on the plus or minus strand!
					 * Problem: Transforming coordinates from negative anchors to plus strand
					 * is not possible without knowing if the read originates from the plus
					 * or the minus strand. This is difficult for reads with few matching anchors
					 * and reads that originate from e.g. an inverted translocation.
					 */
					anchor.onRead = positionOnRead;
					anchor.onRef = part->Scores[k].Location.m_Location;
					if (anchor.isReverse) {
						printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, anchor.onRead, anchor.onRead + readPartLength, part->Scores[k].Location.m_Location + readPartLength, part->Scores[k].Location.m_Location, part->Scores[k].Score.f,
								part->Scores[k].Location.isReverse(),
								DP_TYPE_UNFILTERED, anchor.isUnique ? DP_STATUS_OK : DP_STATUS_LOWSCORE);
					} else {
						printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, anchor.onRead, anchor.onRead + readPartLength, part->Scores[k].Location.m_Location, part->Scores[k].Location.m_Location + readPartLength, part->Scores[k].Score.f,
								part->Scores[k].Location.isReverse(),
								DP_TYPE_UNFILTERED, anchor.isUnique ? DP_STATUS_OK : DP_STATUS_LOWSCORE);
					}

					if (k < 3) {
						double const K = 1 / 3;
						double  const lambda = 1.098612;
						double const n = 256;
						double const m = 3000000000;
						double const Y = part->Scores[k].Score.f;
						double pVal = 1 - exp(-K * n * m * exp(-lambda * Y));
						verbose(0, false, "%f (%f) at %llu, ", part->Scores[k].Score.f, pVal, part->Scores[k].Location.m_Location);
					} else if (k == 3) {
						verbose(0, false, "... (%d)", part->numScores());
					}
				}
				verbose(0, true, "");
			} else {
				verbose(0, true, "no hits found");
				printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, positionOnRead, positionOnRead + readPartLength, 0, 0, 0.0f, 0, DP_TYPE_UNFILTERED, DP_STATUS_NOHIT);
			}
		} else {
			verbose(0, true, "too many hits found");
			printDotPlotLine(group->fullRead->ReadId, group->fullRead->name, positionOnRead, positionOnRead + readPartLength, 0, 0, 0.0f, 0, DP_TYPE_UNFILTERED, DP_STATUS_NOHIT);
		}
	}

	Anchor * anchorsRev = new Anchor[maxAnchorNumber];
	int anchorRevIndex = 0;

	/**
	 * Tree contains all mapped read parts + MQ
	 */
	readCoordsTree = new IntervalTree::IntervalTree<int>(treeIntervals);

	/**
	 * Build HSP from sub read mappings (cLIS)
	 * - Sort sub-reads by read position
	 * - Compute max. number of HSP (based on read length)
	 * - Run cLIS algorithm on reference positions to retrieve HSP
	 * - If isUnique, build Interval from sub-read set and compute regression
	 */
	int nIntervals = 0;
	Interval * * intervals = getIntervalsFromAnchors(nIntervals, anchorsFwd, anchorFwdIndex, anchorsRev, anchorRevIndex, group->fullRead);

	verbose(0, true, "\nIntervals after cLIS:");
	for (int i = 0; i < nIntervals; ++i) {
		verbose(0, "", intervals[i]);
	}
	verbose(0, true, "");

	std::sort(intervals, intervals + nIntervals, sortIntervalsInSegment);

	std::vector<IntervalTree::Interval<Interval *> > intervalList;

	verbose(0, true, "\n\nBuilding segments:\n");
	/**
	 * A segment is a list of intervals that are located in on alignment corridor
	 * Intervals that are contained in others will be deleted. All the others
	 * will be added to a segment
	 */
	int const maxMappedSegementCount = nIntervals + 1;
	MappedSegment * segments = new MappedSegment[maxMappedSegementCount];
	size_t segementsIndex = 0;
	for (int i = 0; i < nIntervals; ++i) {
		Interval * interval = intervals[i];
		bool intervalProcessed = false;
		verbose(0, "Current interval: ", interval);

		for (int j = 0; j < segementsIndex && !intervalProcessed; ++j) {
			verbose(1, true, "Checking segment %d", j);

			for (int k = 0; k < segments[j].length && !intervalProcessed; ++k) {
				Interval * processedInterval = segments[j].list[k];

				verbose(2, "Comparing to interval: ", processedInterval);
				if (isContained(interval, processedInterval)) {
					// Interval is contained in already processed interval and therefore ignored
					intervalProcessed = true;
					verbose(3, true, "Is contained. Deleting interval.");
					delete intervals[i];
					intervals[i] = 0;
					interval = 0;
				} else {
					if (isCompatible(interval, processedInterval)) {
						verbose(3, true, "Is compatible");
						// Interval fits corridor of segment
						if (segments[j].length < segments[j].maxLength) {
							segments[j].list[segments[j].length++] = interval;
							intervalProcessed = true;
							intervalList.push_back(IntervalTree::Interval<Interval *>(interval->onReadStart, interval->onReadStop, interval));
						}
					} else {
						verbose(3, true, "Not contained and not compatible");
					}
				}
			}

		}

		if (!intervalProcessed) {
			// Creating new segment for interval
			verbose(1, true, "Creating new segment %d", segementsIndex);
			if (segementsIndex < maxMappedSegementCount) {
				segments[segementsIndex].list[0] = interval;
				segments[segementsIndex++].length = 1;
				intervalList.push_back(IntervalTree::Interval<Interval *>(interval->onReadStart, interval->onReadStop, interval));
			} else {
				Log.Error("Could not create new segment (%d, %d) for read %s", segementsIndex, maxMappedSegementCount, read->name);
			}
		}

	}

	/**
	 * All intervals are either deleted or added to segments.
	 */
	delete[] intervals;
	intervals = 0;

	// Join segments
	if (pacbioDebug) {
		Log.Message("\nSegments identified: %d", segementsIndex);
		for(int i = 0; i < segementsIndex; ++i) {
			Log.Message("Segment %d contains %d intervals", i, segments[i].length);
		}
	}

	IntervalTree::IntervalTree<Interval *> * intervalsTree = new IntervalTree::IntervalTree<Interval *>(intervalList);

	verbose(0, true, "\n\nMerging segments:\n");
	// Final interval list
	intervals = new Interval *[maxSegmentCount + 1];
	nIntervals = 0;

	Interval * * delIntervals = new Interval *[maxSegmentCount + 1];
	int nDelIntervals = 0;


	/**
	 * Joins segments to full length alignments
	 *  - try to identify all segments that fall into an alignment corridor
	 *  - merge segments that are separated by regions with high error rate
	 *  - merge segments that are separated by small indels (only if both segments are long enough to compensate the score penalty of the deletion)
	 *  - Extend unmerged segments to cover the full read an compensate for missed subread mappings
	 *  - Don't merge segments if separated by a read part that maps to a different place on the genome or is inverted (balanced translocation)
	 */
	verbose(0, true, "Joining segments to intervals:");
	for (int i = 0; i < segementsIndex; ++i) {

		// Sort intervals by position on read
		std::sort(segments[i].list, segments[i].list + segments[i].length, sortIntervalsInSegment);

		if (pacbioDebug) {
			Log.Message("Segment %d:", i);
			for(int j = 0; j < segments[i].length; ++j) {
				verbose(0, "Interval: ", segments[i].list[j]);
			}
		}

		Interval * lastInterval = segments[i].list[0];
		extendIntervalStart(lastInterval, 2 * readPartLength);
		bool isFirstInterval = true;
		Interval * currentInterval = 0;

		for (int j = 1; j < segments[i].length; ++j) {
			currentInterval = segments[i].list[j];
			verbose(1, "Last: ", lastInterval);
			verbose(1, "Current: ", currentInterval);
			if (isSameDirection(currentInterval, lastInterval)) {
				verbose(2, true, "Same direction.");
				loc dupLength = 0;
				if (!isDuplication(currentInterval, lastInterval, dupLength)) {

					if(gapOverlapsWithInterval(lastInterval, currentInterval, intervalsTree, read)) {
						/***
						 * Possible translocation
						 */
						verbose(2, true, "Overlap found in other corridor.");
						verbose(2, true, "Saving last interval. Using current to go on.");

						if (isFirstInterval) {
							extendToReadStart(lastInterval, read->length, intervalsTree, read);
							isFirstInterval = false;
						}
						//TODO: close gap on read (imporve)
						extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
						extendIntervalStart(currentInterval, 2 * readPartLength);
						//Add lastInterval
						intervals[nIntervals++] = lastInterval;
						lastInterval = currentInterval;
					} else {
						/**
						 * Insertion, deletion or gap
						 */
						REAL const corridorSize = std::min(4096, std::min(currentInterval->lengthOnRead(), lastInterval->lengthOnRead()));
						verbose(2, true, "IsContained in corridor of %f.", corridorSize);
						if (canSpanDeletionInsertion(currentInterval, lastInterval, corridorSize) && !spansChromosomeBorder(currentInterval, lastInterval)) {
							/**
							 *  Deletion or insertion small enough for alignment without split
							 */
							verbose(2, true, "Not a duplication, not overlapping with other corridor. Merging current and last interval. Deleting current interval.");
							lastInterval = mergeIntervals(lastInterval, currentInterval);
							segments[i].list[j]->isProcessed = true;
							delIntervals[nDelIntervals++] = segments[i].list[j];
						} else {
							/**
							 * Deletion, insertion or gap
							 */
							verbose(2, true, "Diagonal offset between intervals too big (max %f) or spans chromosome border.", corridorSize);
							verbose(2, true, "Saving last interval. Using current to go on.");

							if (isFirstInterval) {
								extendToReadStart(lastInterval, read->length, intervalsTree, read);
								isFirstInterval = false;
							}
							closeGapOnRead(lastInterval, currentInterval, read->length);
							extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
							extendIntervalStart(currentInterval, 2 * readPartLength);
							//Add lastInterval
							intervals[nIntervals++] = lastInterval;
							lastInterval = currentInterval;
						}
					}
				} else {
					/**
					 * Duplication
					 */
					verbose(2, true, "Covers possible duplication.");
					verbose(2, true, "Saving last interval. Using current to go on.");
					if (isFirstInterval) {
						extendToReadStart(lastInterval, read->length, intervalsTree, read);
						isFirstInterval = false;
					}
					closeGapOnRead(lastInterval, currentInterval, read->length);
					// Extension necessary since border can be wrong even if there is no gap on the read
					int maxExtend = std::min(std::max(currentInterval->onReadStart - lastInterval->onReadStop + (int)dupLength, 0), 2 * readPartLength);
					verbose(2, true, "Duplication stats: %d - %d + %d, %d => %d", currentInterval->onReadStart, lastInterval->onReadStop, (int)dupLength, 2 * readPartLength, maxExtend);
//					int maxExtend = 2 * readPartLength;
					extendIntervalStop(lastInterval, maxExtend, read->length);
					extendIntervalStart(currentInterval, maxExtend);

					//Add lastInterval
					intervals[nIntervals++] = lastInterval;
					lastInterval = currentInterval;
				}
			} else {
				/**
				 * Inversion
				 */
				verbose(2, true, "Not in same direction.");
				verbose(2, true, "Saving last interval. Using current to go on.");
				if (isFirstInterval) {
					extendToReadStart(lastInterval, read->length, intervalsTree, read);
					isFirstInterval = false;
				}
				closeGapOnRead(lastInterval, currentInterval, read->length);
				extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
				extendIntervalStart(currentInterval, 2 * readPartLength);
				//Add lastInterval
				intervals[nIntervals++] = lastInterval;
				lastInterval = currentInterval;
			}
		}
		if (isFirstInterval) {
			extendToReadStart(lastInterval, read->length, intervalsTree, read);
			isFirstInterval = false;
		}
		extendIntervalStop(lastInterval, 2 * readPartLength, read->length);
		verbose(2, true, "Extending last interval to read end", intervalsTree);
		extendToReadStop(lastInterval, read->length, intervalsTree, read);
		verbose(1, "Last: ", lastInterval);
		verbose(2, true, "Saving last interval.");
		intervals[nIntervals++] = lastInterval;
	}

	delete[] segments;
	segments = 0;

	delete intervalsTree;
	intervalsTree = 0;

	for (int i = 0; i < nDelIntervals; ++i) {
		delete delIntervals[i];
		delIntervals[i] = 0;
	}
	delete[] delIntervals;
	delIntervals = 0;
	nDelIntervals = 0;

	verbose(0, true, "Joined intervals:");
	for (int i = 0; i < nIntervals; ++i) {
		verbose(1, "", intervals[i]);
	}
	verbose(0, true, "Sorting intervals by read start");
	std::sort(intervals, intervals + nIntervals, sortIntervalsInSegment);


	if (nIntervals > 0) {
		Interval * lastInterval = intervals[0];
		for (int i = 1; i < nIntervals; ++i) {
			Interval * currentIntervl = intervals[i];

			if (currentIntervl->anchorLength > 1) {
				verbose(1, "a:", lastInterval);
				verbose(1, "b: ", currentIntervl);
				if (!isCompatible(lastInterval, currentIntervl) && getDistanceOnRead(lastInterval, currentIntervl) > 0 && (currentIntervl->anchorLength > 2 || lastInterval->anchorLength > 2)) {
					verbose(1, true, "Closing gap between:");
					closeGapOnRead(lastInterval, currentIntervl, read->length);
				} else {
					verbose(1, true, "Skip");
				}
			}

			if (currentIntervl->anchorLength > 1 || lastInterval->anchorLength == 1) {
				lastInterval = currentIntervl;
			}
		}
	}

	/**
	 * Sort intervals by score
	 * Important because, we trust intervals with higher
	 * score more. Thus we align them first. Aligned intervals
	 * are considered fixed. Therefore, all unangliend intervals will
	 * be trimmed in order to not overlap with fixed intervals.
	 */
	verbose(0, true, "Sorting intervals by score");
	std::sort(intervals, intervals + nIntervals, sortIntervalsByScore);

	int readBpCovered = 0;
	verbose(0, true, "\nFinal intervals:");
	for (int i = 0; i < nIntervals; ++i) {
		verbose(1, "", intervals[i]);
		printDotPlotLine(read->ReadId, read->name, intervals[i]->onReadStart, intervals[i]->onReadStop, intervals[i]->onRefStart, intervals[i]->onRefStop, intervals[i]->score, intervals[i]->isReverse,
		DP_TYPE_SEQMENTS_CONS + i, DP_STATUS_OK);
		readBpCovered += intervals[i]->lengthOnRead();
	}

	float aligned = readBpCovered * 1.0f / read->length;
	verbose(0, true, "Intervals cover %.2f%% of read", aligned * 100.0f);
	bool mapped = (Config.getMinResidues() < 1.0) ? aligned > Config.getMinResidues() : readBpCovered > Config.getMinResidues();
	if (!mapped) {
		verbose(0, true, "Clearing intervals -> read unmapped");
		for (int i = 0; i < nIntervals; ++i) {
			if (intervals[i] != 0) {
				delete intervals[i];
				intervals[i] = 0;
			}
		}
		delete[] intervals;
		intervals = 0;
		nIntervals = 0;

	}

	/**
	 * Align final intervals to the reference
	 */
	if (nIntervals != 0) {

		// Since we don't know how many segments of the read we have to align in the
		// end we need temp arrays to store them
		// TODO: remove fixed upper limit
		LocationScore * tmpLocationScores = new LocationScore[nIntervals * 4];
		Align * tmpAlingments = new Align[nIntervals * 4];
		int nTempAlignments = 0;

		read->Calculated = 0;

		Timer tmr;
		if (pacbioDebug) {
			tmr.ST();
		}
		/**
		 * Aligning intervals
		 */
		for (int i = 0; i < nIntervals; ++i) {
			Interval * currentInterval = intervals[i];

			/**
			 * Adjust unanglined intervals to not overlap with already
			 * aligned intervals
			 */
			verbose(0, "Aligning interval: ", currentInterval);
			for (int j = 0; j < nTempAlignments; ++j) {
				Interval * alignedInterval = tmpAlingments[j].mappedInterval;
				verbose(1, "Check overlap with: ", alignedInterval);

				int overlap = getOverlapOnRead(currentInterval, alignedInterval);
				verbose(1, true, "Overlap: %d", overlap);
				if (overlap > 0 && overlap < currentInterval->lengthOnRead() * 0.95f) {
					verbose(0, true, "Adjusting interval");

					if (currentInterval->onReadStart < alignedInterval->onReadStart) {
						shortenIntervalEnd(currentInterval, overlap);
					} else {
						shortenIntervalStart(currentInterval, overlap);
					}

				}
			}
			verbose(0, "New interval: ", currentInterval);

			/**
			 * Convert intervals on reverse strand to forward strand
			 */
			if (currentInterval->onRefStart > currentInterval->onRefStop) {
				loc tmp = currentInterval->onRefStart;
				currentInterval->onRefStart = currentInterval->onRefStop;
				currentInterval->onRefStop = tmp;
			}

			verbose(0, "Aligning interval: ", currentInterval);
			if (!Config.getSkipalign()) {
				alignSingleOrMultipleIntervals(read, currentInterval, tmpLocationScores, tmpAlingments, nTempAlignments);
			} else {
				Log.Message("Skipping alignment computation.");
			}
			if (nTempAlignments > 0) {
				verbose(0, "Aligned interval: ", tmpAlingments[nTempAlignments - 1].mappedInterval);
			}
		}
		if (pacbioDebug) {
			Log.Message("Alignment took %fs", tmr.ET());
		}

		read->AllocScores(tmpLocationScores, nTempAlignments);
		read->Alignments = tmpAlingments;

		delete[] tmpLocationScores;
		tmpLocationScores = 0;

		/**
		 * Process alignments: removed short alignments, choos combination of aligned segments
		 * that have the highest score and cover the read best
		 */
		if (read->Calculated > 0) {
			bool mapped = reconcileRead(group);
			if (mapped) {
				sortRead(group);
			} else {
				verbose(0, true, "%s (%d) not mapped", read, read->length);
			}
			WriteRead(group->fullRead, mapped);
		} else {
			verbose(0, true, "%s (%d) not mapped", read, read->length);
			WriteRead(group->fullRead, false);
		}
	} else {
		verbose(0, true, "%s (%d) not mapped. No candidates found for read: unmapped.", read, read->length);
		//No candidates found for read
		WriteRead(group->fullRead, false);
	}

	delete[] anchorsFwd;
	anchorsFwd = 0;
	delete[] anchorsRev;
	anchorsRev = 0;

	if (intervals != 0) {
		for (int i = 0; i < nIntervals; ++i) {
			if (intervals[i] != 0) {
				delete intervals[i];
				intervals[i] = 0;
			}
		}
		delete[] intervals;
		intervals = 0;
		nIntervals = 0;
	}

	delete readCoordsTree;
	readCoordsTree = 0;
	processTime += overallTmr.ET();

	verbose(0, true, "###########################################");
	verbose(0, true, "###########################################");
	verbose(0, true, "###########################################");
	verbose(0, true, "");
}

void AlignmentBuffer::SaveRead(MappedRead * read, bool mapped) {
	WriteRead(read, mapped);
}

void AlignmentBuffer::WriteRead(MappedRead* read, bool mapped) {

	if (mapped) {
		//Convert mapping position to RefId and position
		for (int i = 0; i < read->Calculated; ++i) {
			//TODO: fix for -n > 1
			//Instead of setting mapped to false set score to 0 and don't print it in the end
			mapped = SequenceProvider.convert(read->Scores[i].Location);
		}
	}

	m_Writer->WriteRead(read, mapped);

	NGM.GetReadProvider()->DisposeRead(read);
}

