/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "ConvexAlignFast.h"

#include <cmath>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <x86intrin.h>

#include "IConfig.h"

//TODO: remove
#define pRef pBuffer1
#define pQry pBuffer2

namespace Convex {

int ConvexAlignFast::NumberOfSetBits(uint32_t i) {
	// Java: use >>> instead of >>
	// C or C++: use uint32_t
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

ConvexAlignFast::ConvexAlignFast(int const stdOutMode,
		float const match,
		float const mismatch,
		float const gapOpen,
		float const gapExtend,
		float const gapExtendMin,
		float const gapDecay):
		defaultMaxBinaryCigarLength(200000), pacbioDebug(false), stdoutPrintAlignCorridor(stdOutMode) {
	mat = match;
	mis = mismatch;
	gap_open_read = gapOpen;
	gap_open_ref = gapOpen;
	gap_ext = gapExtend;
	gap_decay = gapDecay;
	gap_ext_min = gapExtendMin;

	if(pacbioDebug) {
		fprintf(stderr, "mat: %f mis: %f gap_open_read: %f gap_open_ref: %f gap_ext: %f gap_decay: %f gapExtendMin: %f\n", mat, mis, gap_open_read, gap_open_ref, gap_ext, gap_decay, gapExtendMin);
	}

	matrix = new AlignmentMatrixFast(Config.getMaxMatrixSizeMB());

	maxBinaryCigarLength = defaultMaxBinaryCigarLength;
	binaryCigar = new int[maxBinaryCigarLength];

	alignmentId = 0;
}

ConvexAlignFast::~ConvexAlignFast() {

	delete matrix;
	matrix = 0;

	delete[] binaryCigar;
	binaryCigar = 0;
}

int ConvexAlignFast::printCigarElement(char const op, int const length,
		char * cigar, int & cigarOpCount, int cigarMaxLength) {
	int offset = 0;
	offset = sprintf(cigar, "%d%c", length, op);

	cigarOpCount += 1;

	return offset;
}

void ConvexAlignFast::addPosition(Align & result, int & nmIndex, int posInRef, int posInRead,
		int Yi) {
	if (posInRead > 16 && posInRef > 16) {
		if(nmIndex >= result.nmPerPostionLength) {
			fprintf(stderr, "Debug: PositionNM reallocated.\n");

			int const tmpLength = result.nmPerPostionLength * 2;
			PositionNM * tmp = new PositionNM[tmpLength];

			if(result.nmPerPosition != 0) {
				memcpy(tmp, result.nmPerPosition, result.nmPerPostionLength * sizeof(PositionNM));
				delete[] result.nmPerPosition;
				result.nmPerPosition = 0;
			}
			result.nmPerPosition = tmp;
			result.nmPerPostionLength = tmpLength;
		}
		result.nmPerPosition[nmIndex].readPosition = posInRead - 16;
		result.nmPerPosition[nmIndex].refPosition = posInRef - 16;
		result.nmPerPosition[nmIndex].nm = Yi;
		nmIndex += 1;
	}
}

void ConvexAlignFast::checkMdBufferLength(int md_offset, Align& result, int const minDiff) {
	if (md_offset > (int) ((result.maxMdBufferLength * 0.9f)) || minDiff > (result.maxMdBufferLength - md_offset)) {
		result.maxMdBufferLength *= 2;
		char* tmp = new char[result.maxMdBufferLength];
		memcpy(tmp, result.pQry, sizeof(char) * (result.maxMdBufferLength / 2));
		delete[] result.pQry;
		result.pQry = tmp;
		tmp = 0;
//		fprintf(stderr, "Reallocating MD buffer (%d)\n", result.maxMdBufferLength);
	}
}

int ConvexAlignFast::convertCigar(char const * const refSeq, int const refSeqLength, Align & result,
		FwdResults & fwdResults, int const externalQStart,
		int const externalQEnd) {

	//*********************//
	//Inversion detection init
	//*********************//
	uint buffer = 0;
	int posInRef = 0;
	int posInRead = 0;

	//*********************//
	//General init
	//*********************//

	// Number of CIGAR operations. Currently BAM only supports < 64k CIGAR operations
	int cigarOpCount = 0;

	int nmIndex = 0;
	int exactAlignmentLength = 0;

	int finalCigarLength = 0;

	int cigar_offset = 0;
	int md_offset = 0;

	int binaryCigarIndex = fwdResults.alignment_offset;

	result.svType = 0;
	//*********************//
	// Set QStart
	//*********************//
	result.QStart = ((binaryCigar[binaryCigarIndex] >> 4) + externalQStart);
	if (result.QStart > 0) {
		cigar_offset += printCigarElement('S', result.QStart,
				result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);
		finalCigarLength += result.QStart;
	}
	posInRead = binaryCigar[binaryCigarIndex] >> 4;

	//Positions in read and ref for start of alignment
	result.firstPosition.refPosition = posInRef;
	result.firstPosition.readPosition = posInRead; //QStart of aligned sequence, but not for full read (like result.QStart)

	//*********************//
	// Translate CIGAR to char and compute MD
	//*********************//
	int matches = 0;
	int alignmentLength = 0;

	int cigar_m_length = 0;
	int md_eq_length = 0;
	int ref_index = 0;

	int overallMatchCount = 0;

	//uint const maxIndelLength = 5;
	uint const maxIndelLength = 1;

	//Inversion detection Arndt
	int Yi = 0;

	for (int j = binaryCigarIndex + 1; j < (maxBinaryCigarLength - 1); ++j) {
		int cigarOp = binaryCigar[j] & 15;
		int cigarOpLength = binaryCigar[j] >> 4;

		alignmentLength += cigarOpLength;

		switch (cigarOp) {
		case CIGAR_X:
			cigar_m_length += cigarOpLength;

			//Produces: 	[0-9]+(([A-Z]+|\^[A-Z]+)[0-9]+)*
			//instead of: 	[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
			for (int k = 0; k < cigarOpLength; ++k) {
				checkMdBufferLength(md_offset, result, 100);
				md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
				md_eq_length = 0;
				md_offset += sprintf(result.pQry + md_offset, "%c",
						refSeq[ref_index++]);

				buffer = buffer << 1;
				buffer = buffer | 1;

				//				Yi = std::max(0, Yi + 1);
				Yi = NumberOfSetBits(buffer);
				addPosition(result, nmIndex, posInRef++, posInRead++, Yi);
			}

			exactAlignmentLength += cigarOpLength;

			break;
		case CIGAR_EQ:
			cigar_m_length += cigarOpLength;
			md_eq_length += cigarOpLength;
			matches += cigarOpLength;

			overallMatchCount += cigarOpLength;

			for (int k = 0; k < cigarOpLength; ++k) {
				buffer = buffer << 1;
				//				Yi = std::max(0, Yi - 1);
				Yi = NumberOfSetBits(buffer);
				addPosition(result, nmIndex, posInRef++, posInRead++, Yi);
			}
			ref_index += cigarOpLength;

			exactAlignmentLength += cigarOpLength;

			break;
		case CIGAR_D:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);
				finalCigarLength += cigar_m_length;
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('D', cigarOpLength,
					result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);

			checkMdBufferLength(md_offset, result, 100 + cigarOpLength);
			md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
			md_eq_length = 0;
			result.pQry[md_offset++] = '^';

			for (int k = 0; k < cigarOpLength; ++k) {
				result.pQry[md_offset++] = refSeq[ref_index++];
				buffer = buffer << 1;
				if (k < maxIndelLength) {
					buffer = buffer | 1;
					Yi = std::max(0, Yi + 1);
				}
				addPosition(result, nmIndex, posInRef++, posInRead, Yi);
			}

			exactAlignmentLength += cigarOpLength;

			break;
		case CIGAR_I:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);
				finalCigarLength += cigar_m_length;
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('I', cigarOpLength,
					result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);
			finalCigarLength += cigarOpLength;

			for (int k = 0; k < cigarOpLength; ++k) {
				buffer = buffer << 1;
				if (k < maxIndelLength) {
					buffer = buffer | 1;
					Yi = std::max(0, Yi + 1);
				}
				//				addPosition(result, nmIndex, posInRef++, posInRead, Yi);
				posInRead += 1;
			}

			exactAlignmentLength += cigarOpLength;

			break;
		default:
			fprintf(stderr, "Invalid cigar string: %d\n", cigarOp);
			throw 1;
		}
	}

	//*********************//
	//Print last element
	//*********************//
	checkMdBufferLength(md_offset, result, 100);
	md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
	if (cigar_m_length > 0) {
		cigar_offset += printCigarElement('M', cigar_m_length,
				result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);
		finalCigarLength += cigar_m_length;
		cigar_m_length = 0;
	}

	//*********************//
	//Set QEnd
	//*********************//
	result.QEnd = ((binaryCigar[maxBinaryCigarLength - 1] >> 4) + externalQEnd);
	if (result.QEnd > 0) {
		cigar_offset += printCigarElement('S', result.QEnd,
				result.pRef + cigar_offset, cigarOpCount, result.maxBufferLength);
	}
	finalCigarLength += result.QEnd;

//	fprintf(stderr, "CIGAR: %d, MD: %d, Length: %d\n", cigar_offset, md_offset, result.maxBufferLength);
	if(cigar_offset > result.maxBufferLength || md_offset > result.maxMdBufferLength) {
		fprintf(stderr, "CIGAR/MD buffer not long enough (%d %d > %d %d). Please report this!\n", cigar_offset, md_offset, result.maxBufferLength, result.maxMdBufferLength);
		fprintf(stderr, "CIGAR; %s\n", result.pBuffer1);
		fprintf(stderr, "MD; %s\n", result.pBuffer2);
		throw 1;
	}
	result.Identity = matches * 1.0f / alignmentLength;
	result.pRef[cigar_offset] = '\0';
	result.pQry[md_offset] = '\0';

	result.NM = alignmentLength - matches;
	result.alignmentLength = exactAlignmentLength;

//	if (nmPerPositionLength < exactAlignmentLength) {
//		fprintf(stderr, "Alignmentlength (%d) < exactAlingmentlength (%d)\n",
//				nmPerPositionLength, exactAlignmentLength);
//		throw 1;
//	}
	//	fprintf(stderr, "\n==== Matches: %d of %d ====\n", overallMatchCount,
	//			posInRead);

	//Positions in read and ref for end of alignment
	result.lastPosition.refPosition = posInRef;
	result.lastPosition.readPosition = posInRead;
	result.cigarOpCount = cigarOpCount;
	//	extData[edIndex++] = posInRef;
	//	extData[edIndex++] = posInRead; //QEnd of aligned sequence, but not for full read (like result.QEnd)

	return finalCigarLength;

}

bool ConvexAlignFast::revBacktrack(char const * const refSeq,
		char const * const qrySeq, FwdResults & fwdResults, int readId) {

	if (fwdResults.best_read_index <= 0) {
		return false;
	}

	bool validAlignment = true;

	int binaryCigarIndex = maxBinaryCigarLength - 1;

	int cigar_element = CIGAR_S;
	int cigar_element_length = fwdResults.qend;

	int cigarStringLength = fwdResults.qend;

	int x = fwdResults.best_ref_index;
	int y = fwdResults.best_read_index;

//	AlignmentMatrixFast::MatrixElement currentElement;

	char currentElement = 0;

//	while ((currentElement = matrix->getElement(x, y)->direction) != CIGAR_STOP) {
	while ((currentElement = *matrix->getDirection(x, y)) != CIGAR_STOP) {

//		if (x < 0 || y < 0 || x > matrix->getWidth()
//				|| y > matrix->getHeight()) {
//			fprintf(stderr, "Error in backtracking. x, y indexing error.\n");
//			return false;
//		}

		//TODO: add corridor check. backtracking path too close to corridor end
		if (!matrix->validPath(x, y)) {
			if (pacbioDebug) {
				fprintf(stderr, "Corridor probably too small\n");
			}
			return false;
		}

		if (stdoutPrintAlignCorridor == 6) {
			printf("%d\t%d\t%d\t%d\t%d\n", readId, alignmentId, x, y, 2);
		}

		if (currentElement == CIGAR_X || currentElement == CIGAR_EQ) {
			y -= 1;
			x -= 1;

			cigarStringLength += 1;
		} else if (currentElement == CIGAR_I) {
			y -= 1;
			cigarStringLength += 1;
		} else if (currentElement == CIGAR_D) {
			x -= 1;
		} else {
			fprintf(stderr,
					"Error in backtracking. Invalid CIGAR operation found\n");
			return false;
		}

		if (currentElement == cigar_element) {
			cigar_element_length += 1;
		} else {
			binaryCigar[binaryCigarIndex--] = (cigar_element_length << 4
					| cigar_element);

			cigar_element = currentElement;
			cigar_element_length = 1;
		}

		if (binaryCigarIndex < 0) {
			fprintf(stderr, "Error in backtracking. CIGAR buffer not long enough. Please report this!\n");
			throw 1;
		}

	}

	// Add last element to binary cigar
	binaryCigar[binaryCigarIndex--] =
			(cigar_element_length << 4 | cigar_element);

	binaryCigar[binaryCigarIndex--] = ((y + 1) << 4 | CIGAR_S);
	cigarStringLength += (y + 1);

	fwdResults.ref_position = x + 1;
	fwdResults.qstart = (y + 1);
	//qend was set by "forward" kernel
	fwdResults.alignment_offset = binaryCigarIndex + 1;

	if (matrix->getHeight() != cigarStringLength) {
		fprintf(stderr, "Error read length != cigar length: %d vs %d\n",
				matrix->getHeight(), cigarStringLength);
		validAlignment = false;
	}

	return validAlignment;

}

int ConvexAlignFast::GetScoreBatchSize() const {
	return 0;
}
int ConvexAlignFast::GetAlignBatchSize() const {
	return 0;
}

int ConvexAlignFast::BatchAlign(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, Align * const results,
		void * extData) {

	throw "Not implemented";

	fprintf(stderr, "Unsupported alignment mode %i\n", mode);
	return 0;
}

int ConvexAlignFast::SingleAlign(int const mode, CorridorLine * corridorLines, int const corridorHeight, char const * const refSeq, char const * const qrySeq, Align & align, int const externalQStart, int const externalQEnd, void * extData) {

	alignmentId = align.svType;

	align.svType = 0;
	align.Score = -1.0f;

	int finalCigarLength = -1;

//	try {

	int const refLen = strlen(refSeq);
	int const qryLen = strlen(qrySeq);

	bool allocated = matrix->prepare(refLen, qryLen, corridorLines, corridorHeight);

	if (allocated) {
		align.pBuffer2[0] = '\0';

		FwdResults fwdResults;

		// Debug: rscript convex-align-vis.r
		if (stdoutPrintAlignCorridor == 6) {
			printf("%d\t%d\t%d\t%d\t%d\n", mode, alignmentId, refLen, qryLen, -1);
		}

		AlignmentMatrixFast::Score score = fwdFillMatrixSSESimple(refSeq, qrySeq, fwdResults, mode);

		if (maxBinaryCigarLength < qryLen) {
			maxBinaryCigarLength = qryLen + 1;
			delete[] binaryCigar;
			binaryCigar = 0;
			binaryCigar = new int[maxBinaryCigarLength];
		}

		bool validAlignment = revBacktrack(refSeq, qrySeq, fwdResults, mode);
		if (validAlignment) {
			finalCigarLength = convertCigar(refSeq + fwdResults.ref_position, refLen - fwdResults.ref_position, align, fwdResults, externalQStart, externalQEnd);
			align.PositionOffset = fwdResults.ref_position;
			align.Score = score;

			/**
			 * Check if Qstart clipping was caused by N in reference
			 */
			int const ntestLength = 100;
			int nCount = 0;
			int probeCount = 0;
			for (int k = fwdResults.ref_position; k > std::max(0, fwdResults.ref_position - ntestLength); --k) {
				/**
				 * Decode from SequenceProvider decodes N to X
				 */
				if (refSeq[k] == 'X') {
					nCount += 1;
				}
				probeCount += 1;
			}

			if (nCount > (probeCount * 0.8f)) {
				align.setBitFlag(0x1);
			}
			/**
			 * Check if QEnd clipping was caused by N in reference
			 */
			nCount = 0;
			probeCount = 0;
			for (int k = align.lastPosition.refPosition; k < std::min(align.lastPosition.refPosition + ntestLength, refLen - fwdResults.ref_position); ++k) {
				/**
				 * Decode from SequenceProvider decodes N to X
				 */
				if (refSeq[fwdResults.ref_position + k] == 'X') {
					nCount += 1;
				}
				probeCount += 1;
			}
			if (nCount > (probeCount * 0.8f)) {
				align.setBitFlag(0x1);
			}

		} else {
//			matrix->printMatrix(refSeq, qrySeq);
//			fprintf(stderr, "%d, %d, %d, %d\n", fwdResults.best_read_index, fwdResults.best_ref_index, fwdResults.ref_position, fwdResults.alignment_offset);
			if (pacbioDebug) {
				fprintf(stderr, "Could not backtrack alignment with score %f\n", score);
			}

			align.Score = -1.0f;
			finalCigarLength = -1;
		}
		if (stdoutPrintAlignCorridor == 6) {
			printf("%d\t%d\t%d\t%d\t%d\n", mode, alignmentId, (int) score, finalCigarLength, -3);
		}
	}
//	} catch (...) {
//		fprintf(stderr, "Exception in singlealign\n");
//		align.Score = -1.0f;
//		finalCigarLength = -1;
//	}
	matrix->clean();

	if (maxBinaryCigarLength != defaultMaxBinaryCigarLength) {
		maxBinaryCigarLength = defaultMaxBinaryCigarLength;
		delete[] binaryCigar;
		binaryCigar = 0;
		binaryCigar = new int[maxBinaryCigarLength];
	}

	return finalCigarLength;
}

int ConvexAlignFast::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & align,
		void * extData) {

	fprintf(stderr, "SingleAlign not implemented");
	throw "Not implemented";
}

//#define DEBUG_SSE

#ifdef DEBUG_SSE
#define SZ 100000
float** scorematrix=0;
float** upcellmatrix=0;
float** diagcellmatrix=0;
float** leftcellmatrix=0;
float** maxcellmatrix=0;
int** dirmatrix=0;
int** runmatrix=0;
int** insrunmatrix=0;
int** delrunmatrix=0;

int** alloci()
{
	int** ptr=new int*[SZ];
	for(int i=0;i<SZ;++i)
	{
		ptr[i]=new int[SZ];
	}
	return ptr;	
}

float** allocf()
{
	float** ptr=new float*[SZ];
	for(int i=0;i<SZ;++i)
	{
		ptr[i]=new float[SZ];
	}
	return ptr;	
}

#endif


AlignmentMatrixFast::Score ConvexAlignFast::fwdFillMatrix(char const * const refSeq,
		char const * const qrySeq, FwdResults & fwdResult, int readId) {

	#ifdef DEBUG_SSE
	if(scorematrix==0)
	{
		scorematrix=allocf();
		upcellmatrix=allocf();
		diagcellmatrix=allocf();
		leftcellmatrix=allocf();
		maxcellmatrix=allocf();
		dirmatrix=alloci();
		runmatrix=alloci();
		insrunmatrix=alloci();
		delrunmatrix=alloci();
		fprintf(stderr,"Allocated debug matrices %dx%d\n",SZ,SZ);
	}
	#endif


	AlignmentMatrixFast::Score curr_max = -1.0f;

	for (int y = 0; y < matrix->getHeight(); ++y) {

		matrix->prepareLine(y);

		int xOffset = matrix->getCorridorOffset(y);

		// Debug: rscript convex-align-vis.r
		if (stdoutPrintAlignCorridor == 6) {
			printf("%d\t%d\t%d\t%d\t%d\n", readId, alignmentId, xOffset, y, 0);
			printf("%d\t%d\t%d\t%d\t%d\n", readId, alignmentId,
					xOffset + matrix->getCorridorLength(y), y, 1);
		}

		char const read_char_cache = qrySeq[y];

		for (int x = xOffset; x < (xOffset + matrix->getCorridorLength(y));
				++x) {
			if (x >= matrix->getWidth() || x < 0) {
				continue;
			}

//			AlignmentMatrixFast::MatrixElement const & diag =
//					matrix->getElementConst(x - 1, y - 1);
			AlignmentMatrixFast::Score diag_score = matrix->getScore(x - 1, y - 1);
			AlignmentMatrixFast::MatrixElement const & up = *matrix->getElement(x,
					y - 1);
			AlignmentMatrixFast::MatrixElement const & left = *matrix->getElement(
					x - 1, y);

			bool const eq = read_char_cache == refSeq[x];
			AlignmentMatrixFast::Score const diag_cell = diag_score
					+ ((eq) ? mat : mis);

			AlignmentMatrixFast::Score up_cell = 0;
			AlignmentMatrixFast::Score left_cell = 0;

			int ins_run = 0;
			int del_run = 0;

			if (up.direction == CIGAR_I) {
				ins_run = up.indelRun;
				if (up.score == 0) {
					up_cell = 0;
				} else {
					up_cell = up.score
							+ std::min(gap_ext_min,
									gap_ext + ins_run * gap_decay);
				}
			} else {
				up_cell = up.score + gap_open_read;
			}

			if (left.direction == CIGAR_D) {
				del_run = left.indelRun;
				if (left.score == 0) {
					left_cell = 0;
				} else {
					left_cell = left.score
							+ std::min(gap_ext_min,
									gap_ext + del_run * gap_decay);
				}
			} else {
				left_cell = left.score + gap_open_ref;
			}

			//find max
			AlignmentMatrixFast::Score max_cell = 0;
			max_cell = std::max(left_cell, max_cell);
			max_cell = std::max(diag_cell, max_cell);
			max_cell = std::max(up_cell, max_cell);

			AlignmentMatrixFast::MatrixElement * current = matrix->getElementEdit(x,
					y);
			char & currentDirection = *matrix->getDirection(x, y);

			if (del_run > 0 && max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = del_run + 1;
			} else if (ins_run > 0 && max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = ins_run + 1;
			} else if (max_cell == diag_cell) {
				current->score = max_cell;
				if (eq) {
					current->direction = CIGAR_EQ;
					currentDirection = CIGAR_EQ;
				} else {
					current->direction = CIGAR_X;
					currentDirection = CIGAR_X;
				}
				current->indelRun = 0;
			} else if (max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = 1;
			} else if (max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = 1;
			} else {
				current->score = 0;
				current->direction = CIGAR_STOP;
				currentDirection = CIGAR_STOP;
				current->indelRun = 0;
			}


			#ifdef DEBUG_SSE
			if(x<SZ && y<SZ)
			{
				scorematrix[x][y]=current->score;
				dirmatrix[x][y]=current->direction;
				runmatrix[x][y]=current->indelRun;
				upcellmatrix[x][y]=up_cell;
				leftcellmatrix[x][y]=left_cell;
				diagcellmatrix[x][y]=diag_cell;
				maxcellmatrix[x][y]=max_cell;
				insrunmatrix[x][y]=ins_run;
				delrunmatrix[x][y]=del_run;
			}
			#endif

			//printf("x=%d, y=%d, score=%f, dir=%d, run=%d\n",x,y,current->score,current->direction,current->indelRun);

			if (max_cell > curr_max) {
				curr_max = max_cell;
				fwdResult.best_ref_index = x;
				fwdResult.best_read_index = y;
				fwdResult.max_score = curr_max;
			}

		}

	}
	fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
	if (matrix->getHeight() == 0) {
		fwdResult.best_read_index = fwdResult.best_ref_index = 0;
	}

	return curr_max;
}

AlignmentMatrixFast::Score ConvexAlignFast::FastfwdFillMatrix(char const * const refSeq,
		char const * const qrySeq, FwdResults & fwdResult, int readId) {

	AlignmentMatrixFast::Score curr_max = -1.0f;

	for (int y = 0; y < matrix->getHeight(); ++y) {

		matrix->prepareLine(y);

		int xOffset = matrix->getCorridorOffset(y);

		char const read_char_cache = qrySeq[y];

		int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

		for (int x = std::max(0, xOffset); x < xMax; ++x) {

			AlignmentMatrixFast::Score diag_score = matrix->getElementUp(x - 1, y - 1)->score;
			AlignmentMatrixFast::MatrixElement const & up = *matrix->getElementUp(x,y - 1);
			AlignmentMatrixFast::MatrixElement const & left = *matrix->getElementCurr(x - 1, y);

			bool const eq = read_char_cache == refSeq[x];
			AlignmentMatrixFast::Score const diag_cell = diag_score
					+ ((eq) ? mat : mis);

			AlignmentMatrixFast::Score up_cell = 0;
			AlignmentMatrixFast::Score left_cell = 0;

			int ins_run = 0;
			int del_run = 0;

			if (up.direction == CIGAR_I) {
				ins_run = up.indelRun;
				if (up.score == 0) {
					up_cell = 0;
				} else {
					up_cell = up.score
							+ std::min(gap_ext_min,
									gap_ext + ins_run * gap_decay);
				}
			} else {
				up_cell = up.score + gap_open_read;
			}

			if (left.direction == CIGAR_D) {
				del_run = left.indelRun;
				if (left.score == 0) {
					left_cell = 0;
				} else {
					left_cell = left.score
							+ std::min(gap_ext_min,
									gap_ext + del_run * gap_decay);
				}
			} else {
				left_cell = left.score + gap_open_ref;
			}

			//find max
			AlignmentMatrixFast::Score max_cell = 0;
			max_cell = std::max(left_cell, max_cell);
			max_cell = std::max(diag_cell, max_cell);
			max_cell = std::max(up_cell, max_cell);

			AlignmentMatrixFast::MatrixElement * current = matrix->getElementEditCurr(x,y);

			char & currentDirection = *matrix->getDirection(x, y);

			if (del_run > 0 && max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = del_run + 1;
			} else if (ins_run > 0 && max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = ins_run + 1;
			} else if (max_cell == diag_cell) {
				current->score = max_cell;
				if (eq) {
					current->direction = CIGAR_EQ;
					currentDirection = CIGAR_EQ;
				} else {
					current->direction = CIGAR_X;
					currentDirection = CIGAR_X;
				}
				current->indelRun = 0;
			} else if (max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = 1;
			} else if (max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = 1;
			} else {
				current->score = 0;
				current->direction = CIGAR_STOP;
				currentDirection = CIGAR_STOP;
				current->indelRun = 0;
			}

			if (max_cell > curr_max) {
				curr_max = max_cell;
				fwdResult.best_ref_index = x;
				fwdResult.best_read_index = y;
				fwdResult.max_score = curr_max;
			}

		}

	}
	fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
	if (matrix->getHeight() == 0) {
		fwdResult.best_read_index = fwdResult.best_ref_index = 0;
	}

	return curr_max;
}


#define SIMD_LEVEL 4

typedef union
{
	__m128 v;
	float a[SIMD_LEVEL];
} SSEFloat;

float ssat(__m128 d, int i)
{
	SSEFloat f;
	f.v=d;
	return f.a[i];
}

AlignmentMatrixFast::Score ConvexAlignFast::fwdFillMatrixSSESimple(char const * const refSeq,
		char const * const qrySeq, FwdResults & fwdResult, int readId) {

	#ifdef DEBUG_SSE
	fwdFillMatrix(refSeq,qrySeq,fwdResult,readId);
	#endif

	AlignmentMatrixFast::Score curr_max = -1.0f;

	__m128 mat_sse = _mm_set1_ps(mat);
	__m128 mis_sse = _mm_set1_ps(mis);
	__m128 CIGAR_I_SSE = _mm_set1_ps(CIGAR_I);
	__m128 CIGAR_D_SSE = _mm_set1_ps(CIGAR_D);
	__m128 CIGAR_STOP_SSE = _mm_set1_ps(CIGAR_STOP);
	__m128 CIGAR_X_SSE = _mm_set1_ps(CIGAR_X);
	__m128 CIGAR_EQ_SSE = _mm_set1_ps(CIGAR_EQ);

	__m128 zero_sse = _mm_set1_ps(0);
	__m128 one_sse = _mm_set1_ps(1);

	__m128 gap_ext_min_sse = _mm_set1_ps(gap_ext_min);
	__m128 gap_ext_sse = _mm_set1_ps(gap_ext);
	__m128 gap_decay_sse = _mm_set1_ps(gap_decay);
	__m128 gap_open_read_sse = _mm_set1_ps(gap_open_read);

	for (int y = 0; y < matrix->getHeight(); ++y) {

		matrix->prepareLine(y);

		int xOffset = matrix->getCorridorOffset(y);

		char const read_char_cache = qrySeq[y];
		__m128 read_char_cache_sse = _mm_set1_ps((float)read_char_cache);

		int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

		for (int x = std::max(0, xOffset); x < xMax-SIMD_LEVEL; x+=SIMD_LEVEL) {


			//COMPUTE DIAGONAL SCORE
			__m128 ref_char_cache_sse = _mm_setr_ps((float)refSeq[x],(float)refSeq[x+1],(float)refSeq[x+2],(float)refSeq[x+3]);
			__m128 eq_sse = _mm_cmpeq_ps(read_char_cache_sse,ref_char_cache_sse);


			__m128 diag_score_sse;
			__m128 diag_cell_sse;
			__m128 up_score_sse;
			__m128 up_dir_sse;
			__m128 up_run_sse;


			//COMPUTE UP CELL
			if(x<1 || y==0 || x-1 < matrix->lastCorridor.offset || x+3 >= (matrix->lastCorridor.offset + matrix->lastCorridor.length) )
			{
				const AlignmentMatrixFast::MatrixElement& up0 = *matrix->getElementUp(x,     y - 1);
				const AlignmentMatrixFast::MatrixElement& up1 = *matrix->getElementUp(x + 1, y - 1);
				const AlignmentMatrixFast::MatrixElement& up2 = *matrix->getElementUp(x + 2, y - 1);
				const AlignmentMatrixFast::MatrixElement& up3 = *matrix->getElementUp(x + 3, y - 1);
			
				diag_score_sse = _mm_setr_ps(matrix->getElementUp(x - 1, y - 1)->score,up0.score,up1.score,up2.score);
				diag_cell_sse = _mm_add_ps(diag_score_sse,
                                                          _mm_or_ps(_mm_and_ps(eq_sse, mat_sse),_mm_andnot_ps(eq_sse, mis_sse)));

				up_score_sse = _mm_setr_ps(up0.score,up1.score,up2.score,up3.score);
				up_dir_sse =  _mm_setr_ps(up0.direction,up1.direction,up2.direction,up3.direction);
				up_run_sse =  _mm_setr_ps(up0.indelRun,up1.indelRun,up2.indelRun,up3.indelRun);
			} else {
				const AlignmentMatrixFast::MatrixElement& up0 = *matrix->getElementUpUnprotected(x,     y - 1);
				const AlignmentMatrixFast::MatrixElement& up1 = *matrix->getElementUpUnprotected(x + 1, y - 1);
				const AlignmentMatrixFast::MatrixElement& up2 = *matrix->getElementUpUnprotected(x + 2, y - 1);
				const AlignmentMatrixFast::MatrixElement& up3 = *matrix->getElementUpUnprotected(x + 3, y - 1);
			
				diag_score_sse = _mm_setr_ps(matrix->getElementUpUnprotected(x - 1, y - 1)->score,up0.score,up1.score,up2.score);
				diag_cell_sse = _mm_add_ps(diag_score_sse,
                                                          _mm_or_ps(_mm_and_ps(eq_sse, mat_sse),_mm_andnot_ps(eq_sse, mis_sse)));

				up_score_sse = _mm_setr_ps(up0.score,up1.score,up2.score,up3.score);
				up_dir_sse =  _mm_setr_ps(up0.direction,up1.direction,up2.direction,up3.direction);
				up_run_sse =  _mm_setr_ps(up0.indelRun,up1.indelRun,up2.indelRun,up3.indelRun);
			}	

			//if (up.direction == CIGAR_I) {
			//	ins_run = up.indelRun;
			//	if (up.score == 0) {
			//		up_cell = 0;
			//	} else {
			//		up_cell = up.score
			//				+ std::min(gap_ext_min,
			//						gap_ext + ins_run * gap_decay);
			//	}
			//} else {
			//	up_cell = up.score + gap_open_read;
	                //}
	
			__m128 up_dir_eq_cigari_sse=_mm_cmpeq_ps(up_dir_sse,CIGAR_I_SSE);
			__m128 up_score_eq_zero_sse=_mm_cmpeq_ps(up_score_sse,zero_sse);
			__m128 up_cell_sse=_mm_or_ps( _mm_andnot_ps(up_dir_eq_cigari_sse,_mm_add_ps(up_score_sse,gap_open_read_sse) ),
                                            _mm_and_ps(up_dir_eq_cigari_sse,
                                                      _mm_andnot_ps(up_score_eq_zero_sse,
                                                                 _mm_add_ps(up_score_sse, 
                                                                            _mm_min_ps(gap_ext_min_sse,
                                                                                       _mm_add_ps(gap_ext_sse,
                                                                                                   _mm_mul_ps(up_run_sse,gap_decay_sse) ) )  )  )  ) );			


			//__m128 left_dir_eq_cigard_sse=_mm_cmpeq_ps(left_dir_sse,CIGAR_I_SSE);


			//AlignmentMatrixFast::Score left_cell = 0;

			//if (left.direction == CIGAR_D) {
			//	del_run = left.indelRun;
			//	if (left.score == 0) {
			//		left_cell = 0;
			//	} else {
			//		left_cell = left.score
			//				+ std::min(gap_ext_min,
			//						gap_ext + del_run * gap_decay);
			//	}
			//} else {
			//	left_cell = left.score + gap_open_ref;
			//}

			__m128 curr_dir_sse=CIGAR_STOP_SSE;
			__m128 curr_run_sse=zero_sse;
	
			__m128 max_cell_sse=_mm_max_ps(zero_sse,_mm_max_ps(up_cell_sse,diag_cell_sse));


			//} else if (max_cell == up_cell) {
			//	current->score = max_cell;
			//	current->direction = CIGAR_I;
			//	currentDirection = CIGAR_I;
			//	current->indelRun = 1;

			__m128 cmp_c_sse=_mm_cmpeq_ps(max_cell_sse,up_cell_sse);
			
			curr_dir_sse=_mm_or_ps( _mm_andnot_ps(cmp_c_sse,curr_dir_sse),
                                                _mm_and_ps(cmp_c_sse,CIGAR_I_SSE ) );

			curr_run_sse=_mm_or_ps( _mm_andnot_ps(cmp_c_sse,curr_run_sse),
                                                _mm_and_ps(cmp_c_sse,one_sse) );


			//} else if (max_cell == diag_cell) {
			//	current->score = max_cell;
			//	if (eq) {
			//		current->direction = CIGAR_EQ;
			//		currentDirection = CIGAR_EQ;
			//	} else {
			//		current->direction = CIGAR_X;
			//		currentDirection = CIGAR_X;
			//	}
			//	current->indelRun = 0;
			__m128 cmp_b_sse=_mm_cmpeq_ps(max_cell_sse,diag_cell_sse);
			
			curr_dir_sse=_mm_or_ps( _mm_andnot_ps(cmp_b_sse,curr_dir_sse),
                                                _mm_and_ps(cmp_b_sse,
                                                           _mm_or_ps( _mm_and_ps(eq_sse,CIGAR_EQ_SSE),
                                                                      _mm_andnot_ps(eq_sse,CIGAR_X_SSE) ) ) );

			curr_run_sse=_mm_or_ps( _mm_andnot_ps(cmp_b_sse,curr_run_sse),
                                                _mm_and_ps(cmp_b_sse,zero_sse) );


			//} else if (ins_run > 0 && max_cell == up_cell) {
			//	current->score = max_cell;
			//	current->direction = CIGAR_I;
			//	currentDirection = CIGAR_I;
			//	current->indelRun = ins_run + 1;
			__m128 cmp_a_sse=_mm_and_ps(_mm_cmpgt_ps(up_run_sse,zero_sse), 
                                                    _mm_cmpeq_ps(max_cell_sse,up_cell_sse));

			curr_dir_sse=_mm_or_ps( _mm_andnot_ps(cmp_a_sse,curr_dir_sse),
                                                _mm_and_ps(cmp_a_sse, CIGAR_I_SSE) );

			curr_run_sse=_mm_or_ps( _mm_andnot_ps(cmp_a_sse,curr_run_sse),
                                                _mm_and_ps(cmp_a_sse, _mm_add_ps(up_run_sse,one_sse) ) );


			SSEFloat score_t;
			score_t.v=max_cell_sse;

			SSEFloat dir_t;
			dir_t.v=curr_dir_sse;

			SSEFloat run_t;
			run_t.v=curr_run_sse;

			AlignmentMatrixFast::MatrixElement left = *matrix->getElementCurr(x-1,y);
			for(int j=0;j<SIMD_LEVEL;++j)
			{
				float left_cell=0;
				if (left.direction == CIGAR_D) {
					if (left.score == 0) {
						left_cell = 0;
					} else {
						left_cell = left.score
								+ std::min(gap_ext_min,
										gap_ext + left.indelRun * gap_decay);
					}
				} else {
					left_cell = left.score + gap_open_ref;
				}

				AlignmentMatrixFast::MatrixElement& current = *matrix->getElementCurrUnprotected(x+j,y);
				char & currentDirection = *matrix->getDirection(x+j, y);

				current.score = score_t.a[j];
				current.direction = dir_t.a[j];
				current.indelRun = run_t.a[j];
				currentDirection = dir_t.a[j];

				if(left_cell>=current.score)
				{
					if(left.indelRun > 0)
					{
						current.score=left_cell;
						currentDirection=CIGAR_D;
						current.direction=CIGAR_D;
						current.indelRun=left.indelRun+1;
					} else if ( left_cell>current.score || currentDirection==CIGAR_STOP || (currentDirection==CIGAR_I && ssat(up_run_sse,j)<=0)  ) {
						current.score=left_cell;
						currentDirection=CIGAR_D;
						current.direction=CIGAR_D;
						current.indelRun=1;
					}
				}

				#ifdef DEBUG_SSE
				if(x+j<SZ && y<SZ)
				{
					if(current.score!=scorematrix[x+j][y] || current.direction!=dirmatrix[x+j][y] || current.indelRun!=runmatrix[x+j][y])
					{
						printf("x=%d, y=%d\n",x+j,y);
						printf("score:   %f normal, %f sse\n",scorematrix[x+j][y],current.score);
						printf("dir:     %d normal, %d sse\n",dirmatrix[x+j][y],current.direction);
						printf("run:     %d normal, %d sse\n",current.indelRun,current.indelRun);
						printf("upcell   %f normal, %f sse\n",upcellmatrix[x+j][y],ssat(up_cell_sse,j));
						printf("leftcell %f normal, %f sse\n",leftcellmatrix[x+j][y],left_cell);
						printf("diagcell %f normal, %f sse\n",diagcellmatrix[x+j][y],ssat(diag_cell_sse,j));
						printf("maxcell  %f normal, %f sse\n",maxcellmatrix[x+j][y],ssat(max_cell_sse,j));
						printf("ldelrun  %d normal, %d sse\n",delrunmatrix[x+j][y],left.indelRun);
						printf("uinsrun  %d normal, %f sse\n",insrunmatrix[x+j][y],ssat(up_run_sse,j));
						printf("\n");
						printf("CIGAR_EQ=%d CIGAR_X=%d CIGAR_I=%d CIGAR_D=%d CIGAR_STOP=%d\n",CIGAR_EQ,CIGAR_X,CIGAR_I,CIGAR_D,CIGAR_STOP);
						while(true){}
					}
				}
				#endif

				if (current.score > curr_max) {
					curr_max = current.score;
					fwdResult.best_ref_index = x+j;
					fwdResult.best_read_index = y;
					fwdResult.max_score = curr_max;
				}


				left=current;
			}

		}


		for (int x = std::max(std::max(0,xOffset),xMax-SIMD_LEVEL-8); x < xMax; ++x) {
		//for (int x = std::max(0,xOffset); x < xMax; ++x) {

			AlignmentMatrixFast::Score diag_score = matrix->getElementUp(x - 1, y - 1)->score;
			AlignmentMatrixFast::MatrixElement const & up = *matrix->getElementUp(x,y - 1);
			AlignmentMatrixFast::MatrixElement const & left = *matrix->getElementCurr(x - 1, y);

			bool const eq = read_char_cache == refSeq[x];
			AlignmentMatrixFast::Score const diag_cell = diag_score
					+ ((eq) ? mat : mis);

			AlignmentMatrixFast::Score up_cell = 0;
			AlignmentMatrixFast::Score left_cell = 0;

			int ins_run = 0;
			int del_run = 0;

			if (up.direction == CIGAR_I) {
				ins_run = up.indelRun;
				if (up.score == 0) {
					up_cell = 0;
				} else {
					up_cell = up.score
							+ std::min(gap_ext_min,
									gap_ext + ins_run * gap_decay);
				}
			} else {
				up_cell = up.score + gap_open_read;
			}

			if (left.direction == CIGAR_D) {
				del_run = left.indelRun;
				if (left.score == 0) {
					left_cell = 0;
				} else {
					left_cell = left.score
							+ std::min(gap_ext_min,
									gap_ext + del_run * gap_decay);
				}
			} else {
				left_cell = left.score + gap_open_ref;
			}

			//find max
			AlignmentMatrixFast::Score max_cell = 0;
			max_cell = std::max(left_cell, max_cell);
			max_cell = std::max(diag_cell, max_cell);
			max_cell = std::max(up_cell, max_cell);

			AlignmentMatrixFast::MatrixElement * current = matrix->getElementEditCurr(x,y);
			char & currentDirection = *matrix->getDirection(x, y);


			if (del_run > 0 && max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = del_run + 1;
			} else if (ins_run > 0 && max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = ins_run + 1;
			} else if (max_cell == diag_cell) {
				current->score = max_cell;
				if (eq) {
					current->direction = CIGAR_EQ;
					currentDirection = CIGAR_EQ;
				} else {
					current->direction = CIGAR_X;
					currentDirection = CIGAR_X;
				}
				current->indelRun = 0;
			} else if (max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = 1;
			} else if (max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = 1;
			} else {
				current->score = 0;
				current->direction = CIGAR_STOP;
				currentDirection = CIGAR_STOP;
				current->indelRun = 0;
			}


			if (max_cell > curr_max) {
				curr_max = max_cell;
				fwdResult.best_ref_index = x;
				fwdResult.best_read_index = y;
				fwdResult.max_score = curr_max;
			}

		}

	}

	fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
	if (matrix->getHeight() == 0) {
		fwdResult.best_read_index = fwdResult.best_ref_index = 0;
	}

	return curr_max;
}

/*
AlignmentMatrixFast::Score ConvexAlignFast::FastUnrolledfwdFillMatrix(char const * const refSeq,
		char const * const qrySeq, FwdResults & fwdResult, int readId) {

	AlignmentMatrixFast::Score curr_max = -1.0f;

	for (int y = 0; y < matrix->getHeight(); ++y) {

		matrix->prepareLine(y);
		int xOffset = matrix->getCorridorOffset(y);
		char const read_char_cache = qrySeq[y];
		int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

		for (int x = std::max(0, xOffset); x < xMax; ++x) {
			AlignmentMatrixFast::Score diag_score;
			AlignmentMatrixFast::MatrixElement const * p_up;
			AlignmentMatrixFast::MatrixElement const * p_left;
			AlignmentMatrixFast::MatrixElement * current;

			if(y==0 ||
			  (x-1 <= (matrix->getCorridorOffset(y-1)) ) ||
			  (x >= (matrix->getCorridorOffset(y-1)+matrix->getCorridorLength(y-1)) )  ) {
				diag_score = matrix->getElementUp(x - 1, y - 1)->score;
				p_up = matrix->getElementUp(x,y - 1);
				p_left = matrix->getElementCurr(x - 1, y);
				current = matrix->getElementEditCurr(x,y);
			} else {
				diag_score = matrix->getElementUpFast(x - 1, y - 1)->score;
				p_up = matrix->getElementUpFast(x,y - 1);
				p_left = matrix->getElementCurrFast(x - 1, y);
				current = matrix->getElementEditCurrFast(x,y);
			}
			AlignmentMatrixFast::MatrixElement const & up=*p_up;
			AlignmentMatrixFast::MatrixElement const & left=*p_left;
			char & currentDirection = *matrix->getDirectionCurrFast(x, y);

			bool const eq = read_char_cache == refSeq[x];
			AlignmentMatrixFast::Score const diag_cell = diag_score
					+ ((eq) ? mat : mis);

			AlignmentMatrixFast::Score up_cell = 0;
			AlignmentMatrixFast::Score left_cell = 0;

			int ins_run = 0;
			int del_run = 0;

			if (up.direction == CIGAR_I) {
				ins_run = up.indelRun;
				if (up.score == 0) {
					up_cell = 0;
				} else {
					up_cell = up.score
							+ std::min(gap_ext_min,
									gap_ext + ins_run * gap_decay);
				}
			} else {
				up_cell = up.score + gap_open_read;
			}

			if (left.direction == CIGAR_D) {
				del_run = left.indelRun;
				if (left.score == 0) {
					left_cell = 0;
				} else {
					left_cell = left.score
							+ std::min(gap_ext_min,
									gap_ext + del_run * gap_decay);
				}
			} else {
				left_cell = left.score + gap_open_ref;
			}

			//find max
			AlignmentMatrixFast::Score max_cell = 0;
			max_cell = std::max(left_cell, max_cell);
			max_cell = std::max(diag_cell, max_cell);
			max_cell = std::max(up_cell, max_cell);

			if (del_run > 0 && max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = del_run + 1;
			} else if (ins_run > 0 && max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = ins_run + 1;
			} else if (max_cell == diag_cell) {
				current->score = max_cell;
				if (eq) {
					current->direction = CIGAR_EQ;
					currentDirection = CIGAR_EQ;
				} else {
					current->direction = CIGAR_X;
					currentDirection = CIGAR_X;
				}
				current->indelRun = 0;
			} else if (max_cell == left_cell) {
				current->score = max_cell;
				current->direction = CIGAR_D;
				currentDirection = CIGAR_D;
				current->indelRun = 1;
			} else if (max_cell == up_cell) {
				current->score = max_cell;
				current->direction = CIGAR_I;
				currentDirection = CIGAR_I;
				current->indelRun = 1;
			} else {
				current->score = 0;
				current->direction = CIGAR_STOP;
				currentDirection = CIGAR_STOP;
				current->indelRun = 0;
			}

			if (max_cell > curr_max) {
				curr_max = max_cell;
				fwdResult.best_ref_index = x;
				fwdResult.best_read_index = y;
				fwdResult.max_score = curr_max;
			}

		}

	}
	fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
	if (matrix->getHeight() == 0) {
		fwdResult.best_read_index = fwdResult.best_ref_index = 0;
	}

	return curr_max;
}
*/

int ConvexAlignFast::BatchScore(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results,
		void * extData) {
	throw "Not implemented";
	return 0;

}

}

// g++ -I ../include/ ConvexAlignFast.cpp -o ConvexAlignFast && ./ConvexAlignFast
//int main(int argc, char **argv) {
//
//	IAlignment * aligner = new Convex::ConvexAlignFast(0);
//
////	char const * ref =
////			"TCAAGATGCCATTGTCCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCTCCTGCC";
//	char const * ref =
//			"AGATACATTGTTACTAACTAAGGCCCTACACCAGCAAAGATACATTGTTACTAAGTAAGGCCCTGCACCAACACAG";
//////	char const * qry =
//////			"TCAAGATGCCAATTGTCCCCCGGCCCCCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCC";
////	char const * qry =
////				"CGGGGCCACGGCCACCGCTGCCCTGCC";
////
//	char const * qry =
//			"CCAGAGCAAGTTCCCGGGGGGTGGCGCATAGATGCAGCGGTCTGGGTTCTGTGGCTGGGGGAGTTGGCTTGGGGTCCGTGGCTTGGGTCAGTCTTCAAGAGCCACCCTGGGTGGCACTCCAGGTTCTCTTGGTCTGGGGGGAATTACCCTGACCCCAGCGGTTTCACGGCCTCCTCCCCACTCAGCCTGGGGGAGAGTCCAGGGTCGGCCTGTCCCCCTCCCCCCACTCCTGTCACTATCAGTCCCCTGTGCTCAGCTTACTGGCAGGGTTCCTCCTGGCTGTCCCCTCCTGCCTCCAGCGCTCTTTCCTCAGGTGTCCACGAGCTCGTCCTCATCCCTTTCTGGTCCTGCTCAGATGCCGCCTGAGGCTCCCAAAACAGGGTCTCACTCTGTCACCCAGGGTGGAGTGCAGTGGTGCAATCCAGCTCACTGCATCCTTGACCTCCCAGGCTCAAGCGACCTCTGCCAGCGAGCTTGTTCAATTCCTGCACCAACACGATACATTGTTCTAACTAGGCCCTACACACAGCAAAGATACATTGTACTTAAGTAAGGGCCCTGCACCAACACAGATACATTGTTACTAAGTAGAGGCCCTGCACCAACACACGATACACTTGTTACTAACTAAGGCCCTGCACCAACACAGATACGTTAGTTACTAAGTAAGGCCCTGCACCAACACAGATACGTTGTTACTAAATAAGGCCCTGTACCAACACAGATACATTGTTACTAACTAAGGCCCTGCACCAACACAGATACTTGTTACTAAGTAAGGCCCTGCTACCAACACAGATACATTGTTACTAACTAAGCCCTGCACCAGCCACAGATACATTGTTACTAAGGGCCCTACACCAGCACAGATACATTGTTACTAAGGCCCTGCACCAACACAGATACGTTGTTACTAAGTGAGCCCTGCACCAACACAGATACATTGTTACTAAGTAAGGCCCTGCACCAGCACGAGATACATTGTTACTAAGGCCCTACACCAGCACAGGATAACATTGTCTACTAACTAAGGCCCTGCACGCAACAAAGATACGTATGTTACTAAGTAAGGCCCTGCACCAACACAGATACATTGTTACTAACTAAGGCCCCATGCACCAACACAGATAATTGTTACAAGGCCCTGCACCAACACAGATACATTGTTACTAACTAAGGCCCTGCACCAACACAGATACGTTGTTACTAAGTAAGGCCCTTGCACCAACACAATACATTGTTACTAAATAAGGCCCTGTACCAACATTCAGATACATTGTTACTAACTAAGGCCTGCACCAACACAGATACGTTGTTACTAAGTAAGGCCCTGTACCAACCAGGACTACATTGTTACTAACTAAGGCCCTGCACCACACAGATCGTTGTTTACTAAGTAAGGCCCTGGCACCAACACAGATACATTGTTACTAACTAAGGCCCTGCACAACACAATACATTGTTACTAAGGCCCTGCACCAACACAGATACATTGTTACTAAGTAAGGCCCTGCACCAACACGATACATTGTTACTACTAAGGCCCTGCACCAACACAAATACGTTCTTACTAAGTAAGGCCCTGCACCAACACAGATACATTGTACTAACTAAGGCCCTGTACAACACAGATACGTTGTTACTAAGTAAGGCCCTGCACCAACACAGACTACATTGTTACTAGGCCCTGCACCAAACACAGATACATTGTACTAAGTAAGGCCCTGCAGCCAACTACCAGATGACATTGTTACTAACTAAGGCCCTGCACCAGCACAGATATATTGTTACTAAAGGCCCTGCACCAACCTACAGATACATTGTTACTAACTAAGGCCTGCACCAACAAGATATGTTGTTACTAAGGCCCTGCACTCACACAGATACATTGTTACTAAGTAAGGGCCCTGCACCAACACAGATACATTGTTACTAACTAAAGGCCCCTGCACCAACACAGATAGTTGTTACTAAGTATGGCCCTCACCAACACAGATACATTGTTACTAAGGCCCTGCACCAACACAGATACATGTTACTAAGTAAGGCCCTTGCACCAACACAGATACATTTGTTACTAACTAAGGCCCTGCACCAACACAGATACATTGTACTAAGGCCCTACACCGAGAATGGATACATTGTCACATGGTTACCTAAAGCCTTGCACCAACATGGACACATCATTACTAAACTAAGGCCCTGCACCAACACAGATACATTGTTACTAAGGCCCTGCACCAACACAGATACATTGTTACTAAGATAAGAGCCCTGCACCAACACAGAGTACATTGTTACTAACAAGGCCCTCACCAGCACAGATATATTGTTACTAAGGCCCTGACTCAACACAGATAATTGTTACTAACTCAAGGCCCTGCACCAACCAGATATGTTGTTACTAAGGCCCTGCACCCCAACACAGATACATTGTTACTAAGGCCCGCACCAACACAAGATACATTGTTACTAAGTAAGGCCCTGCACCAACACAGATACGTTGTTACTAACTAAGGCCCTGCACCAACACAGATACGTTGTTACTCGAAGTAAGGCCCTGCACCAATCACAATACATTGTTACTAAGGCCCTGCACCAAACAGATACATTGTTACTAAGTAAGGCCCTGCACCATACACAGATGCATTTTAGCTAACTAAGGCCCTGTACCAACACAGATACATTGTTACTAGCTAAGGCCCTGCACCAACACAAGATACATTGTTACTAAGTAAGGCCCTGTACCAACACAGATATATTGTTACTAACTAAGGCCCTGCACCAACACAGATACGGTTGTTACTAAGTAAGGC";
//	fprintf(stderr, "Ref:  %s\n", ref);
//	fprintf(stderr, "Read: %s\n", qry);
//
//	Align align;
//
//	align.pBuffer1 = new char[strlen(ref) * 4];
//	align.pBuffer2 = new char[strlen(ref) * 4];
//
//	int result = aligner->SingleAlign(0, 20, ref, qry, align, 0);
//	fprintf(stderr, "Alignment terminated with: %d\n", result);
//	fprintf(stderr, "Score: %f\n", align.Score);
//	fprintf(stderr, "Cigar: %s\n", align.pBuffer1);
//	fprintf(stderr, "MD: %s\n", align.pBuffer2);
//
//	delete aligner;
//	aligner = 0;
//
//	return 0;
//}
