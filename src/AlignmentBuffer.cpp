/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
 */

#include "AlignmentBuffer.h"

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <float.h> //Eclipse

#include "Timing.h"
#include "EndToEndAffine.h"
#include "StrippedSW.h"
#include "OutputReadBuffer.h"

using std::vector;

bool AlignmentBuffer::first = true;

void AlignmentBuffer::flush() {
//	DoRun();
//	nReads = 0;
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
//	if (argos) {
//		SaveRead(read, read->hasCandidates());
//	} else {
//		if (!read->hasCandidates() || read->mappingQlty < min_mq) {
//			//If read has no CMRs or mapping quality is lower than min mapping quality, output unmapped read
//			//read->clearScores(-1);
//			SaveRead(read, false);
//		} else {
//			Log.Debug(512, "READ_%d\tALGN_BUFFER\tCMR_%d %f (location %llu) added to alignment buffer at position %d", read->ReadId, scoreID, read->Scores[scoreID].Score.f, read->Scores[scoreID].Location.m_Location, nReads);
//			//add alignment computations to buffer. if buffer is full, submit to CPU/GPU
//			reads[nReads].scoreId = scoreID;
//			reads[nReads++].read = read;
//			if (nReads == batchSize) {
//				DoRun();
//				nReads = 0;
//			}
//		}
//	}
	throw "Not implemented AlignmentBuffer::addRead";
}

void AlignmentBuffer::DoRun() {

//	int count = nReads;
//
//	if (count > 0) {
//		Log.Debug(32, "INFO\tALGN\tSubmitting %d alignment computations.", count);
//		Timer tmr;
//		tmr.ST();
//		for (int i = 0; i < count; ++i) {
//			MappedRead * cur_read = reads[i].read;
//			int scoreID = reads[i].scoreId;
//
//			assert(cur_read->hasCandidates());
//
//			//Initialize
//			if (cur_read->Scores[scoreID].Location.isReverse()) {
//				qryBuffer[i] = cur_read->RevSeq;
//
//				if (cur_read->Paired != 0) {
//					m_DirBuffer[i] = !(cur_read->ReadId & 1);
//				} else {
//					m_DirBuffer[i] = 1;
//				}
//
//			} else {
//				qryBuffer[i] = cur_read->Seq;
//				if (cur_read->Paired != 0) {
//					m_DirBuffer[i] = cur_read->ReadId & 1; //0 if first pair
//				} else {
//					m_DirBuffer[i] = 0;
//				}
//			}
//
//			//decode reference sequence
//			if (!SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
//							cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen)) {
////							cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor)) {
//				Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
//				//Log.Warning("Read sequence: %s", cur_read->Seq);
//				memset(const_cast<char *>(refBuffer[i]), 'N', refMaxLen);
//			}
////			//decode reference sequence
////			SequenceProvider.DecodeRefSequence(const_cast<char *>(refBuffer[i]), 0,
////					cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), refMaxLen);
//
//			//initialize arrays for CIGAR and MD string
//			Log.Error("Should not be here. AlignmentBuffer DoRun");
//			Fatal();
//			static int const qryMaxLen = Config.GetInt("qry_max_len");
//			alignBuffer[i].pBuffer1 = new char[std::max(1, qryMaxLen) * 4];
//			alignBuffer[i].pBuffer2 = new char[std::max(1, qryMaxLen) * 4];
//			*(int*) alignBuffer[i].pBuffer1 = 0x212121;
//			*(int*) alignBuffer[i].pBuffer2 = 0x212121;
//
//			//Log.Message("Ref:  %s\nRead: %s", refBuffer[i], qryBuffer[i]);
//
//		}
//
//		//start alignment
//		int aligned = aligner->BatchAlign(alignmode | (std::max(outputformat, 1) << 8), count, refBuffer, qryBuffer, alignBuffer,
//				(m_EnableBS) ? m_DirBuffer : 0);
//
//		Log.Debug(32, "INFO\tALGN\t%d alignments computed (out of %d)", aligned, count);
//
//		if (aligned != count)
//		Log.Error("Error aligning outputs (%i of %i aligned)", aligned, count);
//
//		//process results
//		for (int i = 0; i < aligned; ++i) {
//			MappedRead * cur_read = reads[i].read;
//			int scoreID = reads[i].scoreId;
//			int id = cur_read->ReadId;
//
//			assert(cur_read->hasCandidates());
//			cur_read->Scores[scoreID].Location.m_Location += alignBuffer[i].PositionOffset - (corridor >> 1);
//
//			cur_read->Alignments[scoreID] = alignBuffer[i];
//
//			Log.Debug(2048, "READ_%d\tALGN_DETAILS\tCMR_%d\t%f\t%f\t%llu\t%.*s\t%s", cur_read->ReadId, scoreID, cur_read->Scores[scoreID].Score.f, alignBuffer[i].Identity, alignBuffer[i].NM, refMaxLen, refBuffer[i], qryBuffer[i]);
//
//			if ((cur_read->Calculated - 1) == scoreID) {
//
//				debugAlgnFinished(cur_read);
//
//				SaveRead(cur_read);
//			}
//
//		}
////		overallTime += tmr.ET();
//	} else {
//		Log.Debug(1, "INFO\tALGN\tEmpty buffer submitted.");
//	}
	throw "Not implemented AlignmentBuffer:DoRun()";
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

CorridorLine * getCorridorEndpointsWithAnchors(Interval const * interval,
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

	return corridorLines;
}

char * extractReferenceSequenceForAlignment(Interval const*& interval,
		int & refSeqLength) {
	char * refSeq = 0;
	refSeqLength = interval->onRefStop - interval->onRefStart + 1;
	if (refSeqLength > 0) {
		//TODO: check why decoded/SequenceProvider writes outside of refSeqLen: This makes + 100 necessary
		refSeq = new char[refSeqLength + 100];
		//decode reference refSeqLength
		if (!SequenceProvider.DecodeRefSequenceExact(refSeq, interval->onRefStart, refSeqLength, 0)) {
			//Log.Warning("Could not decode reference for alignment (read: %s): %llu, %d", cur_read->Scores[scoreID].Location.m_Location - (corridor >> 1), cur_read->length + corridor, cur_read->name);
			Log.Warning("Could not decode reference for alignment");
			delete[] refSeq;
			refSeq = 0;
		}
	}
	return refSeq;
}

Align * AlignmentBuffer::computeAlignment(Interval const * interval,
		int corridor, char * const readSeq, size_t const readLength,
		int const externalQStart, int const externalQEnd, int fullReadLength,
		MappedRead const * const read, bool realign, bool const fullAlignment,
		bool const shortRead = false) {

	if (pacbioDebug) {
		Log.Message("Alignment type: realign %d, full %d, shortRead %d", realign, fullAlignment, shortRead);
	}

	static int alignmentId = 0;

	bool validAlignment = false;

	Align * align = new Align();

	try {

	#ifdef TEST_ALIGNER
		Align * alignFast = new Align();
	#endif

		int refSeqLen = 0;
		char * refSeq = extractReferenceSequenceForAlignment(interval, refSeqLen);

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
			align->pBuffer1 = new char[readLength * 4];
			align->pBuffer2 = new char[readLength * 4];
			align->pBuffer1[0] = '\0';
			align->pBuffer2[0] = '\0';

	#ifdef TEST_ALIGNER
			alignFast->pBuffer1 = new char[readLength * 4];
			alignFast->pBuffer2 = new char[readLength * 4];
			alignFast->pBuffer1[0] = '\0';
			alignFast->pBuffer2[0] = '\0';
	#endif


			int corridorMultiplier = 1;
			while (!validAlignment
					&& (corridor * corridorMultiplier) <= maxCorridorSize
					&& retryCount-- > 0) {

				//Local alignment
				if (pacbioDebug) {
					Log.Message("Aligning %d bp to %d bp (corridor %d)", readLength, refSeqLen, corridor * corridorMultiplier);
					Log.Message("Ref: %s", refSeq);
					Log.Message("Read: %s", readSeq);
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
						corridorLines = getCorridorLinear(corridor * corridorMultiplier, readSeq,
								corridorHeight);
					} else {
						if(corridorMultiplier < 3 && !realign && interval->anchorLength > 0) {
							corridorLines = getCorridorEndpointsWithAnchors(interval,
									corridorMultiplier, refSeq, readSeq, corridorHeight, externalQStart, readPartLength, fullReadLength, realign);
						} else {
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
				int const cigarLength = aligner->SingleAlign(alignmentId, corridorLines,
						corridorHeight, (char const * const ) refSeq,
						(char const * const ) readSeq, *align, externalQStart, externalQEnd, 0);

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
					Log.Message("CIGAR: %s", align->pBuffer1);
					Log.Message("MD:    %s", align->pBuffer2);
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
			Log.Error("Could not align reference sequence for read %s.", read->name);
		}

		if (validAlignment) {
			if (pacbioDebug) {
				Log.Message("%d of %d bp successfully aligned with score %f and identity %f", alignedBp, readLength, align->Score, align->Identity);
			}
			NGM.Stats->alignmentCount += 1;
		} else {
			align->Score = -1.0f;
		}
	} catch (...) {
		Log.Message("Warning: could not compute alignment for read %s", read->name);
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

bool isIntervalInCorridor(REAL k, REAL d, REAL corridor,
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
//	Log.Message("Stop: %llu <= %llu <= %llu", testee.onRefStop, lowerRefStop, upperRefStop);

//TODO: add max distance on read

	return inCorridor;
}

Interval * AlignmentBuffer::toInterval(Anchor const & anchor) {
	if (intervalBufferIndex >= 1000) {
		return 0;
	}

	Interval * interval = new Interval();

	intervalBuffer[intervalBufferIndex++] = interval;

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

bool AlignmentBuffer::isCompatible(Anchor const & anchor,
		Interval const * interval) {
	bool isCompatible = false;

	REAL corridorSize = 512;

	if (anchor.isReverse == interval->isReverse) {
//		isCompatible = (anchor.onRead < interval->onReadStart
//				|| anchor.onRead > interval->onReadStop)
//				&& (anchor.onRef < interval->onRefStart
//						|| anchor.onRef > interval->onRefStop);

		isCompatible = isAnchorInCorridor(interval->m, interval->b,
				corridorSize, anchor);
	}

	return isCompatible;
}

bool AlignmentBuffer::isCompatible(Anchor const & anchor,
		MappedSegment const & segment) {
	bool compatible = false;

	for (int i = 0; i < segment.length && !compatible; ++i) {

		compatible = isCompatible(anchor, segment.list[i]);
	}

	return compatible;
}

// Checks whether a is compatible (is located in the "corridor"
// of b.
bool AlignmentBuffer::isCompatible(Interval const * a, Interval const * b) {

	bool isCompatible = false;
	// Trade off (at the moment):
	// Bigger corridor: higher chance that alignment won't span event (score too low)
	// and a part of the reads is not aligned
	// Smaller corridor: reads are split, SVs are not visible in one alignment
	// TODO: after alignment add check if whole read was aligned. if not,
	// realign unanligned part
	REAL corridorSize = 2048;

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

bool AlignmentBuffer::isContained(Interval const * a, Interval const * b) {
	return a->onReadStart >= b->onReadStart && a->onReadStop <= b->onReadStop

	&& a->onRefStart >= b->onRefStart && a->onRefStop <= b->onRefStop
			&& a->isReverse == b->isReverse;
}

//bool AlignmentBuffer::isContainedOnRead(Interval a, Interval b) {
//	if (abs(a.onReadStop - a.onReadStart) > abs(b.onReadStop - b.onReadStart)) {
//		Interval tmp;
//		tmp = a;
//		a = b;
//		b = tmp;
//	}
//	bool isContained = a.onReadStart >= b.onReadStart
//			&& a.onReadStop <= b.onReadStop;
//
//	return isContained;
//}

Interval * AlignmentBuffer::mergeIntervals(Interval * a, Interval * b) {
//	if (a.onReadStart > b.onReadStart) {
//		a.onReadStart = b.onReadStart;
//		a.onRefStart = b.onRefStart;
//
//	}
//	if (a.onReadStop < b.onReadStop) {
//		a.onReadStop = b.onReadStop;
//		a.onRefStop = b.onRefStop;
//	}

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

bool AlignmentBuffer::constructMappedSegements(MappedSegment * segments,
		Interval * interval, size_t & segmentsIndex) {

	if(pacbioDebug) {
		Log.Message("Adding interval to segments");
	}

	for (int i = 0; i < segmentsIndex; ++i) {
		if(pacbioDebug) {
			Log.Message("\tSegment %d", i);
		}
		Interval * * intervals = segments[i].list;
		int initialLength = segments[i].length;

		for (int j = 0; j < initialLength; ++j) {
			if(pacbioDebug) {
				fprintf(stderr, "\t\t");
				intervals[j]->printOneLine();

			}
			if (isContained(interval, intervals[j])) {
				if (pacbioDebug) {
					Log.Message("Interval is contained in %d. Add and exit.", j);
				}
				//Do not add
				return true;
			} else {
				if (pacbioDebug) {
					Log.Message("Interval not contained in %d. Continue looking.", j);
				}

				if(isCompatible(interval, intervals[j])) {
					if (pacbioDebug) {
						Log.Message("Is compatible. Adding interval to segment and exit");
					}
					if (segments[i].length < 1000) {
						intervals[segments[i].length++] = interval;
						return true;
					} else {
						return false;
					}
				} else {
					if (pacbioDebug) {
						Log.Message("Is NOT compatible. Continue looking");
					}
					//Add
				}
			}
		}
	}

	if (pacbioDebug) {
		Log.Message("Creating new segment %d from interval", segmentsIndex);
	}
	segments[segmentsIndex].list[0] = interval;
	segments[segmentsIndex++].length = 1;

	return true;
}

// Overlap on ref coordinates >= readPartLength
// Overlap on read coordinadtes <= readPartLength
// Difference in overlap must be > 0
// a and b must be on same strand!
bool AlignmentBuffer::isDuplication(Interval const * a, int ttt,
		Interval const * b) {

	if (pacbioDebug) {
		Log.Message("\t\t\tisDuplication:");
		fprintf(stderr, "\t\t\t");
		a->printOneLine();
		fprintf(stderr, "\t\t\t");
		b->printOneLine();
	}

	int overlapOnRead = std::max(0, std::min(a->onReadStop, b->onReadStop) - std::max(a->onReadStart, b->onReadStart));
	loc overlapOnRef = 0;

	if (a->isReverse) {
		// if a and b are on reverse strand, refstart and refstop must be switched
		overlapOnRef = std::max(0ll, std::min(a->onRefStart, b->onRefStart) - std::max(a->onRefStop, b->onRefStop));
	} else {
		overlapOnRef = std::max(0ll, std::min(a->onRefStop, b->onRefStop) - std::max(a->onRefStart, b->onRefStart));
	}

	if (pacbioDebug) {
		Log.Message("\t\t\toverlap on read: %d, overlap on ref: %lld", overlapOnRead, overlapOnRef);
	}

	loc overlapDiff = std::max(0ll, overlapOnRef - overlapOnRead);

	return overlapOnRef >= (int) (1.0f * readPartLength)
			&& overlapOnRead <= (int) (readPartLength * 1.0f) && overlapDiff > 0ll;
}

bool sortIntervalsInSegment(Interval const * a, Interval const * b) {
	return a->onReadStart < b->onReadStart;
}

Interval * * AlignmentBuffer::consolidateSegments(MappedSegment * segments,
		size_t segmentsIndex, int & intervalsIndex) {

	if (pacbioDebug) {
		Log.Message("==+==========================");
		Log.Message("=====consolidateSegments=====");
		Log.Message("===+=========================");
	}

	int maxIntervalCount = 0;
	for (int i = 0; i < segmentsIndex; ++i) {
		maxIntervalCount += segments[i].length;

		if(pacbioDebug) {
			Log.Message("Segment %d consists of %d intervals", i, segments[i].length);
		}

	}

	Interval * * intervals = new Interval * [maxIntervalCount + 1];

	for (int i = 0; i < segmentsIndex; ++i) {


		std::sort(segments[i].list, segments[i].list + segments[i].length,
				sortIntervalsInSegment);

		if (pacbioDebug) {
			Log.Message("Segment %d: ", i);
			for(int j = 0; j < segments[i].length; ++j) {
				fprintf(stderr, "\tInterval %d: ", j);
				segments[i].list[j]->printOneLine();
			}
		}

		Interval * lastInterval = segments[i].list[0];
		Interval * currentInterval = 0;

		for (int j = 1; j < segments[i].length; ++j) {
			currentInterval = segments[i].list[j];
			if(pacbioDebug) {
				Log.Message("\tCurrent interval is: %d", j);
			}
			if (isSameDirection(currentInterval, lastInterval)
					&& !isDuplication(currentInterval, 0, lastInterval)) {
				lastInterval = mergeIntervals(currentInterval, lastInterval);
				if(pacbioDebug) {
					Log.Message("\t\tMerging with last interval");
				}
			} else {
				//If not single assigned index
				// Single mapped subreads are added to segment.
				// If they merge with an interval -> included
				// If not -> discarded
				if(!currentInterval->isAssigned) {
					//Add lastInterval
					if(!lastInterval->isAssigned) {
						if(pacbioDebug) {
							fprintf(stderr, "\t\tAdding interval: ");
							lastInterval->printOneLine();
						}
						intervals[intervalsIndex++] = lastInterval;
					}
					lastInterval = currentInterval;
				}
			}
		}
		if(!lastInterval->isAssigned) {
			if(pacbioDebug) {
				fprintf(stderr, "\tAdding interval %d: ", intervalsIndex);
				lastInterval->printOneLine();
			}
			intervals[intervalsIndex++] = lastInterval;
		}

	}
	if (pacbioDebug) {
		Log.Message("=============================");
		Log.Message("===consolidateSegments end===");
		Log.Message("=============================");
	}

	return intervals;
}

/* Takes all anchors found for the read and transforms it in to intervals.
 * Intervals are continues candidate mapping for a read or parts of the
 * read. Contains repeatedly computing cLIS and reconciling of the found
 * intervals. */
Interval * * AlignmentBuffer::infereCMRsfromAnchors(int & intervalsIndex,
		Anchor * allFwdAnchors, int allFwdAnchorsLength, Anchor * allRevAnchors,
		int allRevAnchorsLength, MappedRead * read) {

	if (pacbioDebug) {
		Log.Message("Finding LIS for read %d", read->ReadId);
	}

	//TODO: remove fixed length
	int const maxMappedSegementCount = maxIntervalNumber;
	MappedSegment * segments = new MappedSegment[maxMappedSegementCount];
	size_t segementsIndex = 0;

	//Sort by position on read. Probably not necessary!!
	std::sort(allFwdAnchors, allFwdAnchors + allFwdAnchorsLength,
			sortAnchorOnRead);

	int const maxcLISRunNumber = maxIntervalNumber;
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

			if(lisLength >= 1) {
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

				REAL * regX = new REAL[std::max(2, allFwdAnchorsLength)];
				REAL * regY = new REAL[std::max(2, allFwdAnchorsLength)];
				int pointNumber = 0;

				//Find intervall that covers all anchors best
				Interval * interval = new Interval();
				// This is terrible. It is a result of the overall terrible code design here.
				// It tracks all created intervals and is only used to delete them
				// before processLongRead terminates
				intervalBuffer[intervalBufferIndex++] = interval;

				interval->anchors = new Anchor[lisLength + 1];
				interval->anchorLength = 0;

				for (int i = 0; i < allFwdAnchorsLength; ++i) {
					if (posInLIS >= 0 && i == lis[posInLIS]) {

						interval->anchors[interval->anchorLength++] = allFwdAnchors[lis[posInLIS]];

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

							printDotPlotLine(read->ReadId, read->name,
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

							printDotPlotLine(read->ReadId, read->name,
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
				if(pointNumber == 1) {
					regX[0] = minOnRef;
					regY[0] = minOnRead;
					regX[1] = maxOnRef;
					regY[1] = maxOnRead;
					pointNumber = 2;
				}
				linreg(pointNumber,regX,regY,&m,&b,&r);
				delete[] regX; regX = 0;
				delete[] regY; regY = 0;
				if(pacbioDebug) {
					Log.Message("Regression: m=%.*f b=%.*f r=%.*f\n", DBL_DIG,m,DBL_DIG,b,DBL_DIG,r);
				}

				interval->isReverse = isReverse;
				interval->score = intervalScore;

				interval->onReadStart = minOnRead;
				interval->onReadStop = maxOnRead;

				interval->onRefStart = minOnRef;
				interval->onRefStop = maxOnRef;

				interval->m = m;
				interval->b = b;
				interval->r = r;

				if(interval->m != 0.0) {
					printDotPlotLine(read->ReadId, read->name,
							interval->m, interval->b, r,
							interval->score,
							interval->isReverse,
							DP_TYPE_SEQMENTS_REG + cLISRunNumber, DP_STATUS_NOCOORDS);

//				printDotPlotLine(id, name,
//						-1 * interval.m, -1 * interval.b, r,
//						interval.score,
//						interval.isReverse,
//						DP_TYPE_SEQMENTS_REG + cLISRunNumber * 2 + 1, DP_STATUS_NOCOORDS);
				}

				printDotPlotLine(read->ReadId, read->name,
						interval->onReadStart,
						interval->onReadStop,
						interval->onRefStart,
						interval->onRefStop,
						interval->score,
						interval->isReverse,
						DP_TYPE_SEQMENTS + cLISRunNumber, DP_STATUS_OK);

				//TODO: check chromosome borders

				if (pacbioDebug) {
					Log.Message("New interval:");
					interval->print();
				}

				if(abs(interval->onReadStop - interval->onReadStart) > 0 &&
						llabs(interval->onRefStart - interval->onRefStop) > 0) {
					if(!constructMappedSegements(segments, interval, segementsIndex)) {
						lisLength = 0;
						delete[] lis;
						Log.Error("Too many intervals found (>1000). Ignoring all others");
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
			}

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

	if (allFwdAnchorsLength > 0) {

		Timer anchorTimer;
		anchorTimer.ST();
		if (pacbioDebug) {
			Log.Message("Unassigned anchors (%d):", allFwdAnchorsLength);
			std::cerr << "Elements: (i:onRead:onref:isReverse) ";
			for (int i = 0; i < allFwdAnchorsLength; ++i) {
				std::cerr << i << ":" << allFwdAnchors[i].onRead << ":"
				<< allFwdAnchors[i].onRef << ":"
				<< allFwdAnchors[i].isReverse << ", ";
			}
			std::cerr << std::endl;
		}

		int assigned = 0;
		for(int i = 0; i < allFwdAnchorsLength; ++i) {
			bool added = false;
			for(int j = 0; j < segementsIndex && !added; ++j) {
				if(isCompatible(allFwdAnchors[i], segments[j])) {
					//Add anchor to segement
					added = true;
					assigned += 1;
					addAnchorAsInterval(allFwdAnchors[i], segments[j]);
				}
				//allFwdAnchors[i]
			}
		}

		if(pacbioDebug) {
			Log.Message("Assigned anchors: %d (took %fs", assigned, anchorTimer.ET());
		}
	}

	intervalsIndex = 0;
	Interval * * intervals = consolidateSegments(segments, segementsIndex, intervalsIndex);
	delete[] segments;
	segments = 0;

	// Extend all intervals to account for unmapped anchors at the end
	for(int i = 0; i < intervalsIndex; ++i) {
		Interval * interval = intervals[i];

		// If possible add 512 bp to start end end
		// min to avoid extending beyond the ends of the read
		int addStart = std::min(1 * readPartLength, interval->onReadStart);
		int addEnd = std::min(1 * readPartLength, read->length - interval->onReadStop);

		// Try to align full read, only possible if one fragment, otherwise might align through SVs
		// Could be extended to > 1 fragment: only extend left most and right most fragement into one direction

		int intervalReadLength = (interval->onReadStop - interval->onReadStart);

		if(intervalReadLength >= (read->length * 0.9f)) {
			addStart = interval->onReadStart - 1;
			addEnd = read->length - interval->onReadStop - 1;
		}

		// Try to match the slope of the interval when extending
		double lengthRatio = (interval->onReadStop - interval->onReadStart) * 1.0f / llabs(interval->onRefStop - interval->onRefStart) * 1.0f;

		if(interval->isReverse) {
			interval->onReadStart -= addStart;
			interval->onRefStart += (int) round(addStart / lengthRatio);

			interval->onReadStop += addEnd;
			interval->onRefStop -= (int) round(addEnd / lengthRatio);
		} else {
			interval->onReadStart -= addStart;
			interval->onRefStart -= (int) round(addStart / lengthRatio);

			interval->onReadStop += addEnd;
			interval->onRefStop += (int) round(addEnd / lengthRatio);

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

//int AlignmentBuffer::alignmentCheckForInversion(int const inversionLength,
//		const int refCheckLength, SequenceLocation inversionCheckLocation,
//		uloc inversionMidpointOnRead, const char* const readName,
//		int inversionNumber, char* fullReadSeq) {
//
//	static bool const noLowQualSplit = !Config.getLowQualitySplit();
//
//	int svType = SV_NONE;
//	int const readCheckLength = 50;
//
//	int const fullReadSeqLength = strlen(fullReadSeq);
//
//	// TODO: don't allocate and delete every time. Can be done once
//	// and reused
//	const int readSeqLenght = readCheckLength * 2 + 10;
//	char* readSeq = new char[readSeqLenght];
//	char* revreadSeq = new char[readSeqLenght];
//	memset(readSeq, '\0', readSeqLenght * sizeof(char));
//	memset(revreadSeq, '\0', readSeqLenght * sizeof(char));
//	const int refSeqLength = inversionLength + 2 * refCheckLength;
//	char* refSeq = new char[refSeqLength + 10];
//	memset(refSeq, '\0', (refSeqLength + 10) * sizeof(char));
//	SequenceProvider.DecodeRefSequence(refSeq, 0, inversionCheckLocation.m_Location, refSeqLength);
//
//	// Extract sequence from read.
//	if (pacbioDebug) {
//		Log.Message("CheckLength: %d, Midpoint: %d, readSeqLength: %d", readCheckLength, inversionMidpointOnRead, readSeqLenght);
//	}
//	if (readCheckLength <= inversionMidpointOnRead
//			&& (inversionMidpointOnRead + readCheckLength)
//					< fullReadSeqLength) {
//		strncpy(readSeq,
//				fullReadSeq + inversionMidpointOnRead - readCheckLength,
//				readCheckLength * 2);
//	}
//
////	if (readCheckLength > inversionMidpointOnRead) {
////		strncpy(readSeq, fullReadSeq,
////				std::min(readCheckLength * 2, fullReadSeqLength));
////	} else {
////		if ((inversionMidpointOnRead + readCheckLength) < readSeqLenght) {
////			strncpy(readSeq,
////					fullReadSeq + inversionMidpointOnRead - readCheckLength,
////					readCheckLength * 2);
////		}
////	}
//	if (strlen(readSeq) > 0) {
//		// Reverse complement
//		computeReverseSeq(readSeq, revreadSeq, readCheckLength * 2);
//		if (pacbioDebug) {
//			Log.Message("Ref:      %s", refSeq);Log.Message("Read:     %s", readSeq);Log.Message(
//					"RevRead:  %s", revreadSeq);}
//		if (printInvCandidateFa) {
//			printf(">%s_%d/1\n%s\n", readName, inversionNumber, refSeq);
//			printf(">%s_%d/2\n%s\n", readName, inversionNumber, revreadSeq);
//		}
//		IAlignment* aligner = new EndToEndAffine();
//		float scoreFwd = 0.0f;
//		float scoreRev = 0.0f;
////		static const float minScore = Config.GetFloat(MATCH_BONUS) * 1.0f
////				* readCheckLength / 4.0f;
//		//Formula from Moritz -> p-val > 0.05 for Match:1, mismacth: -4, gap open: -2, gap extend -1
//		const float minScore = 15.0f;
//		aligner->SingleScore(10, 0, refSeq, readSeq, scoreFwd, 0);
//		aligner->SingleScore(10, 0, refSeq, revreadSeq, scoreRev, 0);
//		delete aligner;
//		aligner = 0;
//		if (pacbioDebug) {
//			SequenceProvider.convert(inversionCheckLocation);
//			int len = 0;
//			Log.Message(
//					"Inversion check at position %s:%llu - fwd: %f, rev: %f, min: %f",
//					SequenceProvider.GetRefName(inversionCheckLocation.getrefId(), len),
//					inversionCheckLocation.m_Location, scoreFwd, scoreRev, minScore);
//		}
//		if (scoreRev > scoreFwd && scoreRev > minScore) {
//			svType = SV_INVERSION;
//		} else if (scoreRev < minScore && scoreFwd < minScore
//				&& !noLowQualSplit) {
//
//			svType = SV_TRANSLOCATION;
//		}
//		//		inversionVerified = scoreRev > scoreFwd || scoreFwd < minScore;
//	} else {
//		if (pacbioDebug) {
//			Log.Message("Skipping alignment. Readpart not long enough.");
//		}
//	}
//	delete[] refSeq;
//	refSeq = 0;
//	delete[] readSeq;
//	readSeq = 0;
//	delete[] revreadSeq;
//	revreadSeq = 0;
//
//	return svType;
//}

int AlignmentBuffer::checkForSV(Align const * const align, Interval const * interval, char * fullReadSeq, uloc inversionMidpointOnRef, uloc inversionMidpointOnRead, int const inversionLength, MappedRead * read) {

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
			IAlignment* aligner = new EndToEndAffine();
//			IAlignment* aligner = new StrippedSW();

			float scoreFwd = 0.0f;
			float scoreRev = 0.0f;
			static const float minScore = Config.getScoreMatch() * 1.0f * readCheckLength / 4.0f;
			//Formula from Moritz -> p-val > 0.05 for Match:1, mismacth: -4, gap open: -2, gap extend -1
//			const float minScore = 15.0f;
			aligner->SingleScore(10, 0, refSeq, readSeq, scoreFwd, 0);
			aligner->SingleScore(10, 0, refSeq, revreadSeq, scoreRev, 0);


//			Align alingtest;
//			alingtest.pBuffer1 = new char[10000];
//			alingtest.pBuffer2 = new char[10000];
//			aligner->SingleAlign(10, 0, refSeq, readSeq, alingtest, nullptr);
//			Log.Message("Fwd: %s", alingtest.pBuffer1);
//
//			aligner->SingleAlign(10, 0, refSeq, revreadSeq, alingtest, nullptr);
//			Log.Message("Rev: %s", alingtest.pBuffer1);
			delete aligner;
			aligner = 0;
			if (pacbioDebug) {
				SequenceProvider.convert(inversionCheckLocation);
				int len = 0;
				Log.Message(
						"Inversion check at position %s:%llu - fwd: %f, rev: %f, min: %f",
						SequenceProvider.GetRefName(inversionCheckLocation.getrefId(), len),
						inversionCheckLocation.m_Location, scoreFwd, scoreRev, minScore);
			}
			if (scoreRev >= scoreFwd && scoreRev > minScore) {
				svType = SV_INVERSION;
			} else if (scoreRev < minScore && scoreFwd < minScore && !noLowQualSplit) {

				svType = SV_TRANSLOCATION;
			}
			//		inversionVerified = scoreRev > scoreFwd || scoreFwd < minScore;
		} else {
			if (pacbioDebug) {
				Log.Message("Skipping alignment. Readpart not long enough.");
			}
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
		Interval const * alignedInterval, char * readPartSeq,
		Interval * leftOfInv, Interval * rightOfInv, MappedRead * read) {

	static bool const enabled = Config.getSmallInversionDetection();
	if (!enabled) {
		return SV_NONE;
	}

	//*********************//
	//Detect inversions
	//*********************//

	uloc bestInversionMidpointOnRead = 0;
	uloc bestInversionMidpointOnRef = 0;

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
	mappingLocation.m_Location = alignedInterval->onRefStart
			+ align->PositionOffset;
	SequenceProvider.convert(mappingLocation);

	if (stdoutErrorProfile) {
		int len = 0;
		for (int i = 0; i < align->alignmentLength; ++i) {
			printf("%s\t%llu\t%d\t%s\n", SequenceProvider.GetRefName(mappingLocation.getrefId(), len), mappingLocation.m_Location + align->nmPerPosition[i].refPosition, align->nmPerPosition[i].nm, read->name);
		}
	}

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

Align * AlignmentBuffer::alignInterval(MappedRead const * const read_,
		Interval const * interval, char * const readSeq,
		size_t const readSeqLen, bool const realign, bool const fullAlignment) {

	Align * align = 0;

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
			QStart = read_->length - interval->onReadStop;

		} else {
			QStart = interval->onReadStart;
			QEnd = read_->length - interval->onReadStop;
		}

		if (pacbioDebug) {
			Log.Message("Computing alignment - Start pos: %d, Length: %d", interval->onReadStart, readSeqLen);
		}
		align = computeAlignment(interval, corridor, readSeq, readSeqLen,
				QStart, QEnd, read_->length, read_, realign, fullAlignment,
				false);
	} else {
		if (pacbioDebug) {
			Log.Warning("Tried to align invalid interval:");
			interval->print();
		}
	}
	alignTime += alignTimer.ET();

	return align;
}

char* const AlignmentBuffer::extractReadSeq(const size_t& readSeqLen,
		Interval const * interval, MappedRead* read,
		bool const revComp = false) {
	if (pacbioDebug) {
		Log.Message("Allocating %d bytes", readSeqLen + 1);
	}

	char * readSeq = new char[readSeqLen + 1];
	if (interval->isReverse) {
		read->computeReverseSeq();
		computeReverseSeq(read->Seq + interval->onReadStart, readSeq, readSeqLen);
		readSeq[readSeqLen] = '\0';
	} else {
		strncpy(readSeq, read->Seq + interval->onReadStart, readSeqLen);
		readSeq[readSeqLen] = '\0';
	}
	if (pacbioDebug) {
		Log.Message("ReadSeqLen: %d", readSeqLen);
	}
	if (pacbioDebug) {
		Log.Message("Start pos: %d, Length: %d", interval->onReadStart, readSeqLen);
	}

	if(revComp) {
		char * tmp = new char[readSeqLen + 1];
		computeReverseSeq(readSeq, tmp, readSeqLen);
		delete[] readSeq;
		readSeq = tmp;
		readSeq[readSeqLen] = '\0';
	}

	return readSeq;
}

int AlignmentBuffer::realign(int const svTypeDetected,
		Interval const * interval, Interval * leftOfInv, Interval * rightOfInv,
		MappedRead * read, Align * tmpAling, int & alignIxndex,
		LocationScore * tmpScore, int mq) {

	if (pacbioDebug) {
		Log.Message("***********************************");
		Log.Message("Realigning read!");
		Log.Message("Interval:");
		interval->print();
	}

	float realignScore = 0.0f;
	int svTypeResult = 0;

	/**********************************************/
	//Align part of read that is left of inversion
	/**********************************************/
	if (pacbioDebug) {
		Log.Message("Left Interval:");
		leftOfInv->print();
	}
	size_t readSeqLen = leftOfInv->onReadStop - leftOfInv->onReadStart;
	char * readSeq = extractReadSeq(readSeqLen, leftOfInv, read);

	Align * alignLeft = alignInterval(read, leftOfInv, readSeq, readSeqLen, true, false);
	LocationScore scoreLeft;
	Align * alignRight = 0;
	LocationScore scoreRight;
	Align * alignInv = 0;
	LocationScore scoreInv;

	delete[] readSeq;
	readSeq = 0;

	// Defines the inverted part of the alignment
	// start and end are computed from the alignments
	// left and right to it
	Interval * inv = new Interval();
	if (alignLeft != 0 && alignLeft->Score > 0.0f) {
		alignLeft->MQ = mq;
		//tmpAling[alignIndex] = alignLeft;

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
		if (pacbioDebug) {
			Log.Message("Right Interval:");
			rightOfInv->print();
		}
		readSeqLen = rightOfInv->onReadStop - rightOfInv->onReadStart;
		readSeq = extractReadSeq(readSeqLen, rightOfInv, read);

		alignRight = alignInterval(read, rightOfInv, readSeq, readSeqLen, true, false);
		delete[] readSeq;
		readSeq = 0;

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
			//int inversionLength = inv->onReadStop - inv->onReadStart;
			int inversionLength = llabs(inv->onRefStop - inv->onRefStart);
			int const minInversionLength = Config.getMinInversionLength();
//			int extLength = (int) round((inv->onReadStop - inv->onReadStart) * 0.2f);
//			inv->onReadStart = std::max(0, inv->onReadStart - extLength);
//			inv->onReadStop = std::min(read->length, inv->onReadStop + extLength);
//			inv->onRefStart -= extLength;
//			inv->onRefStop += extLength;

//			if(svTypeDetected == SV_INVERSION) {
			if (pacbioDebug) {
				Log.Message("Length of inversion is %d", inversionLength);
				Log.Message("Inverted Interval:");
				inv->print();
			}
			if(inversionLength > minInversionLength) {
				readSeqLen = inv->onReadStop - inv->onReadStart;
				readSeq = extractReadSeq(readSeqLen, inv, read);
				alignInv = alignInterval(read, inv, readSeq, readSeqLen, true, true);
				delete[] readSeq;
				readSeq = 0;
				if(pacbioDebug && alignInv != 0) {
					Log.Message("Inversion alignment score: %f", alignInv->Score);
				}

				char * revReadSeq = extractReadSeq(readSeqLen, inv, read, true);
				Align * alignInvRev = alignInterval(read, inv, revReadSeq, readSeqLen, true, true);
				delete[] revReadSeq;
				revReadSeq = 0;
				if(pacbioDebug && alignInvRev != 0) {
					Log.Message("Inversion alignment score: %f", alignInvRev->Score);
				}

				if (alignInv != 0 && alignInv->Score > 0.0f && alignInv->getAlignedReadBp(read->length) > minInversionLength && (alignInvRev == 0 || alignInvRev->Score < alignInv->Score)) {
					alignInv->MQ = mq;
					alignInv->svType = SV_INVERSION;

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
					//int len = 0;
					//printf("%s\t%llu\t%llu\n", SequenceProvider.GetRefName(inversionStartRef.getrefId(), len), inversionStartRef.m_Location, inversionStopRef.m_Location);

					svTypeResult = SV_INVERSION;
				} else {
					if (pacbioDebug) {
						Log.Message("Alignment of inverted part failed or score not high enough. Reporting translocation.");
					}
					svTypeResult = SV_TRANSLOCATION;
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
		tmpAling[alignIxndex] = *alignLeft;
		alignLeft->clearNmPerPosition();
		delete alignLeft;
		alignLeft = 0;
		tmpScore[alignIxndex] = scoreLeft;
		alignIxndex += 1;
		read->Calculated += 1;
		tmpAling[alignIxndex] = *alignRight;
		alignRight->clearNmPerPosition();
		delete alignRight;
		alignRight = 0;
		tmpScore[alignIxndex] = scoreRight;
		alignIxndex += 1;
		read->Calculated += 1;

		if(svTypeResult == SV_INVERSION && alignInv != 0) {
			tmpAling[alignIxndex] = *alignInv;
			alignInv->clearNmPerPosition();
			delete alignInv;
			alignInv = 0;
			tmpScore[alignIxndex] = scoreInv;
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
	static float minResidues = Config.getMinResidues();
	static float minIdentity = Config.getMinIdentity();

	if (minResidues <= 1.0f) {
		minResidues = readLength * minResidues;
	}
	return (align->Score > 0.0f) && (align->Identity >= minIdentity) && ((float) (readLength - align->QStart - align->QEnd) >= minResidues);
}


void AlignmentBuffer::alignSingleOrMultipleIntervals(MappedRead * read,
		Interval const * const interval, LocationScore * tmp, Align * tmpAling,
		int & alignIndex) {

//	if (interval->anchorLength > 1) {

	size_t const readSeqLen = interval->onReadStop - interval->onReadStart;
	char * readPartSeq = extractReadSeq(readSeqLen, interval, read);

	Align * align = alignInterval(read, interval, readPartSeq, readSeqLen,
			false, false);
	if (align != 0) {
		if (align->Score > 0.0f) {
			if (pacbioDebug) {
				SequenceLocation seqLoc2;
				seqLoc2.m_Location = interval->onRefStart
						+ align->PositionOffset;
				SequenceProvider.convert(seqLoc2);
				Log.Message("Mapping Pos (%s): %llu -- %llu", read->name, interval->onRefStart + align->PositionOffset, seqLoc2.m_Location);
			}

			Interval * leftOfInv = new Interval();
			Interval * rightOfInv = new Interval();
			int svType = SV_UNKNOWN;
			bool inversionAligned = false;
			svType = detectMisalignment(align, interval, readPartSeq, leftOfInv,
					rightOfInv, read);

			if (svType != SV_NONE) {
				int mq = computeMappingQuality(*align, read->length);
				int assumedSvType = svType;
				svType = realign(svType, interval, leftOfInv, rightOfInv, read,
						tmpAling, alignIndex, tmp, mq);
				if (svType == SV_NONE && assumedSvType == SV_INVERSION) {
					SequenceLocation loc;
					loc.m_Location = leftOfInv->onRefStop;
					SequenceProvider.convert(loc);
					SequenceLocation loc2;
					loc2.m_Location = rightOfInv->onRefStart;
					SequenceProvider.convert(loc2);
					int len = 0;
					if(pacbioDebug) {
						Log.Message("Potential SV detected at %s:%llu-%llu but not taken into account due to errors. (read: %s)", SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location, loc2.m_Location, read->name);
					}
				}
			} else {
				if (pacbioDebug) {
					Log.Message("No SV detected!");
				}
			}
			delete leftOfInv;
			leftOfInv = 0;
			delete rightOfInv;
			rightOfInv = 0;

			if (svType == SV_NONE) {
				if(pacbioDebug) {
					Log.Message("No SV was detected. Keeping single alignment!");
				}
				/**********************************************/
				// No inversion detected
				/**********************************************/
				if (satisfiesConstraints(align, read->length)) {
					align->MQ = computeMappingQuality(*align, read->length);
					tmpAling[alignIndex] = *align;
					align->clearNmPerPosition();
					delete align;
					align = 0;

					tmp[alignIndex].Location.m_Location = interval->onRefStart
					+ tmpAling[alignIndex].PositionOffset;//- (corridor >> 1); handled in computeAlingment
					tmp[alignIndex].Location.setReverse(interval->isReverse);
					tmp[alignIndex].Score.f = tmpAling[alignIndex].Score;

					read->Calculated += 1;
					alignIndex += 1;
				} else {
					align->clearBuffer();
					align->clearNmPerPosition();
					delete align;
					align = 0;
					if (pacbioDebug) {
						Log.Message("Alignment failed");
					}
				}

				// TODO: realign
//			/**********************************************/
//			// Check if significant portion of read was
//			// not aligned
//			/**********************************************/
//			//TODO: add remapping functionality
//			int const maxClipping = readSeqLen * 0.4f;
//			int const minLenght = readPartLength * 2.0f;
//
//			int missingLength = align->firstPosition.readPosition;
//			if(missingLength > maxClipping
//					&& missingLength > minLenght) {
//				if(pacbioDebug) {
//					SequenceLocation loc;
//					loc.m_Location = tmp[alignIndex - 1].Location.m_Location;
//					SequenceProvider.convert(loc);
//					int len = 0;
//					Log.Message("Read: %s (%s:%llu)", read->name, SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location);
//					Log.Message("Skipped %d bases at beginning of %d bp long read (%f)", missingLength, readSeqLen, (missingLength * 100.0f / readSeqLen));
//				}
//			}
//
//			missingLength = readSeqLen - align->lastPosition.readPosition;
//			if(missingLength > maxClipping && missingLength > minLenght) {
//				if(pacbioDebug) {
//					SequenceLocation loc;
//					loc.m_Location = tmp[alignIndex - 1].Location.m_Location;
//					SequenceProvider.convert(loc);
//					int len = 0;
//					Log.Message("Read: %s (%s:%llu)", read->name, SequenceProvider.GetRefName(loc.getrefId(), len), loc.m_Location);
//					Log.Message("Skipped %d bases at end of %d bp long read (%f)", missingLength, readSeqLen, (missingLength * 100.0f / readSeqLen));
//				}
//			}
			} else {
				align->clearBuffer();
				align->clearNmPerPosition();
				delete align;
				align = 0;
			}
		} else {
			align->clearBuffer();
			align->clearNmPerPosition();
			delete align;
			align = 0;

			if (pacbioDebug) {
				Log.Message("Alignment failed");
			}
		}
	}
	delete[] readPartSeq;
	readPartSeq = 0;

//	} else {
//		if (pacbioDebug) {
//			Log.Message("Skipping, not enough anchor support.");
//		}
//	}
}

int AlignmentBuffer::computeMappingQuality(Align const & alignment,
		int readLength) {
	std::vector<IntervalTree::Interval<int> > results;

	readCoordsTree->findOverlapping(alignment.QStart,
			readLength - alignment.QEnd, results);
	int mqSum = 0;
	int mqCount = 0;
	for (int j = 0; j < results.size(); ++j) {
		mqSum += results[j].value;
		mqCount += 1;
	}
	return (int) (mqSum * 1.0f / mqCount);
}

bool sortMappedSegements(Interval * a,
		Interval * b) {
//	return a.value.onReadStart < b.value.onReadStart;
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

void AlignmentBuffer::reconcileRead(ReadGroup * group) {

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

	for(int i = 0; i < bestSegments.size(); ++i) {
		int index = bestSegments[i];

		if (pacbioDebug) {
			Log.Message("\tFragment: %d", index);
		}
		mappedSegements[index]->isProcessed = true;

		if(mappedSegements[index]->score > topFragmentScore) {
			fragmentWithHighestScore = index;
		}
	}
	// Set the alignment with the higest score as primary
	if(bestSegments.size() > 0) {
		read->Alignments[mappedSegements[fragmentWithHighestScore]->id].primary = true;
	}

	bestSegments.clear();
	float secondBestScore = getBestSegmentCombination(maxLength, mappedSegements, bestSegments);
	if(pacbioDebug) {
		Log.Message("Second best score: %f", secondBestScore);
	}

	// Mark remaining as skipped (will not be written to SAM file)
	for (int i = 0; i < mappedSegements.size(); ++i) {
		if(!mappedSegements[i]->isProcessed) {
			read->Alignments[mappedSegements[i]->id].skip = true;
		}
	}

	// DEBUG code: print final list of mappeg segments
	for (int i = 0; i < read->Calculated; ++i) {
		if (!read->Alignments[mappedSegements[i]->id].skip) {
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


	// Delete all mapped segments
	for (int i = 0; i < mappedSegements.size(); ++i) {
		delete mappedSegements[i];
		mappedSegements[i] = 0;
	}
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

void AlignmentBuffer::processShortRead(MappedRead * read) {
	// Read is shorter then anchor length

	if (pacbioDebug) {
		Log.Message("Read too short for long read cLIS");
		Log.Message("Processing ShortRead: %d - %s (lenght %d)", read->ReadId, read->name, read->length);
	}

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

			if (pacbioDebug) {
				Log.Message("Aligning interval: ");
				interval->print();
			}

			int const shortReadCorridor = readPartLength * 1 + 2 * refExtend;

			char * readPartSeq = extractReadSeq(read->length, interval, read);

//			Align * align = alignInterval(read, interval, read->Seq, read->length, false);
			Align * align = computeAlignment(interval, shortReadCorridor, readPartSeq,
					read->length, 0, 0, read->length, read, false, false, true);

			delete[] readPartSeq;
			readPartSeq = 0;

			if (align != 0 && align->Score > 0.0f) {
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
				if (pacbioDebug) {
					Log.Message("Alignment failed");
				}
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

void AlignmentBuffer::processLongReadLIS(ReadGroup * group) {

//	WriteRead(group->fullRead, false);
//	return;

	Timer overallTmr;
	overallTmr.ST();

	if (group->fullRead->length > readPartLength) {

		// Parts of read that were aligned plus MQ of subalignments
		std::vector<IntervalTree::Interval<int> > treeIntervals;

		if (pacbioDebug) {
			Log.Message("Processing LongReadLIS: %d - %s (lenght %d)", group->fullRead->ReadId, group->fullRead->name, group->fullRead->length);
		}

//		float avgGroupScore = (group->bestScoreSum * 1.0f / group->readsFinished) / readPartLength;
//		float minGroupScore = avgGroupScore * 0.5f;

		// TODO: remove fixed length and don't allocate and delete for
		// every read
		int maxAnchorNumber = 10000;
		Anchor * anchorsFwd = new Anchor[maxAnchorNumber];
		int anchorFwdIndex = 0;

		Anchor * anchorsRev = new Anchor[maxAnchorNumber];
		int anchorRevIndex = 0;

		for (int j = 0; j < group->readNumber; ++j) {
			MappedRead * part = group->reads[j];

			int positionOnRead = j * readPartLength;

			treeIntervals.push_back(
					IntervalTree::Interval<int>(positionOnRead,
							positionOnRead + readPartLength,
							part->mappingQlty));

			// Min score required to consider an anchor for the candidate search
			float minScore = (
					(part->numScores() > 0) ?
							part->Scores[0].Score.f * 0.9 : 0.0f)
					/ readPartLength;
//			minScore = std::max(minScore, minGroupScore);
//		minScore = 0.0f;

			// Get all mapping positions from anchors (non overlapping 512bp parts of reads)
			// and convert them to Anchor objects
			// TODO: remove fixed length of 100
			int maxCandidatesPerReadPart = 100;
			for (int k = 0; k < part->numScores(); ++k) {

				if (stdoutPrintScores) {
					printf("%f\n", part->Scores[k].Score.f);
				}

				// If anchor has to many mapping positions -> ignore
				if (part->numScores() < maxCandidatesPerReadPart) {
					if ((part->Scores[k].Score.f / readPartLength) > minScore) {
						// Anchor is valid and will be used
						if (anchorFwdIndex >= (maxAnchorNumber - 1)) {
							Log.Message("Anchor array too small - reallocating.");
							maxAnchorNumber = maxAnchorNumber * 2;
							Anchor * anchorsTmp = new Anchor[maxAnchorNumber];
							memcpy(anchorsTmp, anchorsFwd, anchorFwdIndex * sizeof(Anchor *));
							delete[] anchorsFwd;
							anchorsFwd = anchorsTmp;
							delete[] anchorsRev;
							anchorsRev = new Anchor[maxAnchorNumber];
						}
						Anchor & anchor = anchorsFwd[anchorFwdIndex++];

						anchor.score = part->Scores[k].Score.f;
						anchor.isReverse = part->Scores[k].Location.isReverse();
						anchor.type = DP_STATUS_OK;

						// It would be best to convert reads or read parts that map to the negative strand
						// Immediately to plus strand.
						// If somebody tells me how to do this properly, I'll happily change this.
						// For now Anchors can be on the plus or minus strand!
						// Problem: Transforming coordinates from negative anchors to plus strand
						// is not possible without knowing if the read originates from the plus
						// or the minus strand. This is difficult for reads with few matching anchors
						// and reads that originate from e.g. an inverted translocation.
						anchor.onRead = positionOnRead;
						anchor.onRef = part->Scores[k].Location.m_Location;
						if (anchor.isReverse) {
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
						if(part->Scores[k].Location.isReverse()) {
							printDotPlotLine(group->fullRead->ReadId,
									group->fullRead->name, positionOnRead, positionOnRead + readPartLength,
									part->Scores[k].Location.m_Location + readPartLength,
									part->Scores[k].Location.m_Location,
									part->Scores[k].Score.f,
									part->Scores[k].Location.isReverse(),
									DP_TYPE_UNFILTERED, DP_STATUS_LOWSCORE);
						} else {
							printDotPlotLine(group->fullRead->ReadId,
									group->fullRead->name, positionOnRead, positionOnRead + readPartLength,
									part->Scores[k].Location.m_Location,
									part->Scores[k].Location.m_Location + readPartLength,
									part->Scores[k].Score.f,
									part->Scores[k].Location.isReverse(),
									DP_TYPE_UNFILTERED, DP_STATUS_LOWSCORE);
						}
					}
				} else {
					// Repetitive
					// Still print for visualization
					if(part->Scores[k].Location.isReverse()) {
						printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
								positionOnRead, positionOnRead + readPartLength,
								part->Scores[k].Location.m_Location + readPartLength,
								part->Scores[k].Location.m_Location,
								part->Scores[k].Score.f,
								part->Scores[k].Location.isReverse(),
								DP_TYPE_UNFILTERED,
								DP_STATUS_REPETITIVE);
					} else {
						printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
								positionOnRead, positionOnRead + readPartLength,
								part->Scores[k].Location.m_Location,
								part->Scores[k].Location.m_Location + readPartLength,
								part->Scores[k].Score.f,
								part->Scores[k].Location.isReverse(),
								DP_TYPE_UNFILTERED,
								DP_STATUS_REPETITIVE);
					}
				}
			}

			if (part->numScores() == 0) {
				//No hits found
				if (pacbioDebug) {
					Log.Message("No hits found for part %d", positionOnRead);
				}
				printDotPlotLine(group->fullRead->ReadId, group->fullRead->name,
						positionOnRead, positionOnRead + readPartLength, 0, 0, 0.0f, 0, DP_TYPE_UNFILTERED, DP_STATUS_NOHIT);
			}
		}

		readCoordsTree = new IntervalTree::IntervalTree<int>(treeIntervals);

		int nIntervals = 0;
		Interval * * intervals = infereCMRsfromAnchors(nIntervals, anchorsFwd,
				anchorFwdIndex, anchorsRev, anchorRevIndex, group->fullRead);

		if (nIntervals != 0) {

			if (pacbioDebug) {
				Log.Message("================Intervalls found================");
				for (int i = 0; i < nIntervals; ++i) {
					Interval * interval = intervals[i];
					interval->printOneLine();
				}
				Log.Message("================++++++++++++++++================");
			}

			//Prepare alignment list in read object
			MappedRead * read = group->fullRead;

			// Since we don't know how many segments of the read we have to align in the
			// end we need temp arrays to store them
			// TODO: remove fixed upper limit
			LocationScore * tmpLocationScores = new LocationScore[nIntervals * 4];
			Align * tmpAlingments = new Align[nIntervals * 4];
			int nTempAlignments = 0;

			read->Calculated = 0;

			if (pacbioDebug) {
				Log.Message("========================");
			}
			Timer tmr;
			if (pacbioDebug) {
				tmr.ST();
			}

			/********************************/
			/********** Align CMRs **********/
			/********************************/
			for (int i = 0; i < nIntervals; ++i) {

				if (pacbioDebug) {
					Log.Message("############################################");
					Log.Message("Aligning interval %d", i);
					intervals[i]->print();

				}
				printDotPlotLine(read->ReadId, read->name,
						intervals[i]->onReadStart,
						intervals[i]->onReadStop,
						intervals[i]->onRefStart,
						intervals[i]->onRefStop,
						intervals[i]->score,
						intervals[i]->isReverse,
						DP_TYPE_SEQMENTS_CONS + i, DP_STATUS_OK);

				if (intervals[i]->onRefStart > intervals[i]->onRefStop) {
					loc tmp = intervals[i]->onRefStart;
					intervals[i]->onRefStart = intervals[i]->onRefStop;
					intervals[i]->onRefStop = tmp;
				}

				alignSingleOrMultipleIntervals(read, intervals[i], tmpLocationScores, tmpAlingments, nTempAlignments);

			}
			if (pacbioDebug) {
				Log.Message("Alignment took %fs", tmr.ET());
			}

			if (pacbioDebug) {
				Log.Message("================Intervalls aligned================");
				int bpAligned = 0;

				for(int i = 0; i < read->Calculated; ++i) {
					Align align = tmpAlingments[i];
					Log.Message("%d - %d", align.QStart, read->length - align.QEnd);
					bpAligned += (read->length - align.QStart - align.QEnd);
				}
				Log.Message("Aligned %d bp of %d bp", bpAligned, read->length);
				Log.Message("================++++++++++++++++++=================");
			}

			read->AllocScores(tmpLocationScores, nTempAlignments);
			read->Alignments = tmpAlingments;

			delete[] tmpLocationScores;
			tmpLocationScores = 0;

			if (pacbioDebug) {
				Log.Message("==========================");
			}
			if (read->Calculated > 0) {

				reconcileRead(group);

				sortRead(group);

				WriteRead(group->fullRead, true);
			} else {
				WriteRead(group->fullRead, false);
			}
		} else {
			//No candidates found for read
			WriteRead(group->fullRead, false);
		}

		delete[] anchorsFwd;
		anchorsFwd = 0;
		delete[] anchorsRev;
		anchorsRev = 0;
		delete[] intervals;
		intervals = 0;

		delete readCoordsTree;
		readCoordsTree = 0;
	} else {
		Log.Message("Should not be here. Short read in long read pipeline.");
	}

	for (int i = 0; i < intervalBufferIndex; ++i) {
		delete intervalBuffer[i];
		intervalBuffer[i] = 0;
	}
	intervalBufferIndex = 0;
	processTime += overallTmr.ET();
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

