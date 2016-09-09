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

#ifndef CONVEXALIGN_FAST_H_
#define CONVEXALIGN_FAST_H_

#include "IAlignment.h"

#include <stdint.h>
#include <algorithm>

#include "Types.h"
#include "AlignmentMatrixFast.h"

namespace Convex {

class ConvexAlignFast: public IAlignment {

public:
	ConvexAlignFast(int gpu_id);

	virtual ~ConvexAlignFast();

	virtual int GetScoreBatchSize() const;

	virtual int GetAlignBatchSize() const;

	/**
	 * Not implemented!
	 */
	virtual int BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, float * const results,
			void * extData);

	/**
	 * Not implemented!
	 */
	virtual int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results,
			void * extData);

	virtual int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData);

	/**
	 * Main function. Compute alignment based on "corridor layout"
	 */
	virtual int SingleAlign(int const mode, CorridorLine * corridor,
			int const corridorHeight, char const * const refSeq,
			char const * const qrySeq, Align & result, int const externalQStart,
			int const externalQEnd, void * extData);

private:

	/**
	 * Information that hase to be passed
	 * from fwd to backward and compute cigar step
	 */
	struct FwdResults {
		int best_ref_index;
		int best_read_index;
		int max_score;
		int qend;
		int qstart;
		int ref_position;
		int alignment_offset;
	};

	/**
	 * Scoring function
	 */
	AlignmentMatrixFast::Score mat;
	AlignmentMatrixFast::Score mis;
	AlignmentMatrixFast::Score gap_open_read;
	AlignmentMatrixFast::Score gap_open_ref;
	AlignmentMatrixFast::Score gap_ext;
	AlignmentMatrixFast::Score gap_ext_min;
	AlignmentMatrixFast::Score gap_decay;

	AlignmentMatrixFast * matrix;

	/** Temporary storage for binary
	 * representation of the CIGAR string
	 * produced by backtracking step.
	 */
	int * binaryCigar;

	/**
	 * Max length of binaryCigar
	 */
	int const maxBinaryCigarLength;

	/**
	 * Debug flag
	 * TODO: convert to preprocessor define
	 */
	bool const pacbioDebug;

	/**
	 * Prints debug information for
	 * aligment corridor visualization
	 * to stdout
	 */
	bool const stdoutPrintAlignCorridor;

	/**
	 * Debug
	 */
	int alignmentId;

	/**
	 * Print CIGAR element to cigar string
	 * Helper function used by convertCigar
	 */
	int printCigarElement(char const op, int const length, char * cigar);

	/**
	 * Forward step of alignment
	 * Fills directionMatrix
	 */
	AlignmentMatrixFast::Score fwdFillMatrix(char const * const refSeq,
				char const * const qrySeq, FwdResults & fwdResults, int readId);

	AlignmentMatrixFast::Score FastfwdFillMatrix(char const * const refSeq,
				char const * const qrySeq, FwdResults & fwdResults, int readId);

	AlignmentMatrixFast::Score FastUnrolledfwdFillMatrix(char const * const refSeq,
				char const * const qrySeq, FwdResults & fwdResults, int readId);

	/*AlignmentMatrixFast::Score FastUnrolledfwdFillMatrixMaster(char const * const refSeq,
				char const * const qrySeq, FwdResults & fwdResults, int readId);

	AlignmentMatrixFast::Score FastUnrolledfwdFillMatrixLine(char const * const refSeq,
				char const * const qrySeq, FwdResults & fwdResults, int readId,
				AlignmentMatrixFast::Score diag_score,
				AlignmentMatrixFast::MatrixElement const & up,
				AlignmentMatrixFast::MatrixElement const & left,
				AlignmentMatrixFast::MatrixElement * current,
				char & currentDirection,
				AlignmentMatrixFast::Score & curr_max,
				int x,
				int y);*/

	inline AlignmentMatrixFast::Score FastUnrolledfwdFillMatrixMaster(char const * const refSeq,
			char const * const qrySeq, FwdResults & fwdResult, int readId) {

		AlignmentMatrixFast::Score curr_max = -1.0f;
		int matrix_height=matrix->getHeight();


		{
			int y=0;

			matrix->prepareLine(y);
			int xOffset = matrix->getCorridorOffset(y);
			int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

			for (int x = std::max(0, xOffset); x < xMax; ++x) {
				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
					matrix->empty.score,
					matrix->empty,
					*matrix->getElementCurr(x - 1, y),
					matrix->getElementEditCurrFast(x,y),
					*matrix->getDirectionCurrFast(x, y),
					curr_max,
					x,
					y,
					qrySeq[0]);
			}
		}


		for (int y = 1; y < matrix_height; ++y) {

			matrix->prepareLine(y);
			int xOffset = matrix->getCorridorOffset(y);

			int xMin=std::max(0, xOffset);
			int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

			int xFastMin=std::max(0,matrix->getCorridorOffset(y-1)+1);
			int xFastMax=std::min(matrix->getCorridorOffset(y-1)+matrix->getCorridorLength(y-1),matrix->getWidth());


			char const read_char_cache = qrySeq[y];

			if(xMin>xFastMin || xFastMax>xMax)
			{
				for(int x=xMin; x<xMax; x++)
				{
					ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
						matrix->getElementUp(x - 1, y - 1)->score,
						*matrix->getElementUp(x,y - 1),
						*matrix->getElementCurr(x - 1, y),
						matrix->getElementEditCurrFast(x,y),
						*matrix->getDirectionCurrFast(x, y),
						curr_max,
						x,
						y,
						read_char_cache);
				}			
			} else {
				for(int x=xMin; x<=xFastMin; x++)
				{
					ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
						matrix->getElementUp(x - 1, y - 1)->score,
						*matrix->getElementUp(x,y - 1),
						*matrix->getElementCurr(x - 1, y),
						matrix->getElementEditCurrFast(x,y),
						*matrix->getDirectionCurrFast(x, y),
						curr_max,
						x,
						y,
						read_char_cache);
				}

				for(int x=xFastMin+1; x<xFastMax; x++)
				{
					ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
						matrix->getElementUpFast(x - 1, y - 1)->score,
						*matrix->getElementUpFast(x,y - 1),
						*matrix->getElementCurrFast(x - 1, y),
						matrix->getElementEditCurrFast(x,y),
						*matrix->getDirectionCurrFast(x, y),
						curr_max,
						x,
						y,
						read_char_cache);
				}


				for(int x=xFastMax; x<xMax; x++)
				{
					ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
						matrix->getElementUp(x - 1, y - 1)->score,
						*matrix->getElementUp(x,y - 1),
						*matrix->getElementCurr(x - 1, y),
						matrix->getElementEditCurrFast(x,y),
						*matrix->getDirectionCurrFast(x, y),
						curr_max,
						x,
						y,
						read_char_cache);		
				}			
			}


	 
			/*for (int x = std::max(0, xOffset); x < xMax; ++x) {
				if(y==0 ||
				  (x-1 <= (matrix->getCorridorOffset(y-1)) ) ||
				  (x >= (matrix->getCorridorOffset(y-1)+matrix->getCorridorLength(y-1)) )  ) {
					ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
						matrix->getElementUp(x - 1, y - 1)->score,
						*matrix->getElementUp(x,y - 1),
						*matrix->getElementCurr(x - 1, y),
						matrix->getElementEditCurrFast(x,y),
						*matrix->getDirectionCurrFast(x, y),
						curr_max,
						x,
						y);
				} else {
					ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
						matrix->getElementUpFast(x - 1, y - 1)->score,
						*matrix->getElementUpFast(x,y - 1),
						*matrix->getElementCurrFast(x - 1, y),
						matrix->getElementEditCurrFast(x,y),
						*matrix->getDirectionCurrFast(x, y),
						curr_max,
						x,
						y);
				}
			}*/

		}
		fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
		if (matrix->getHeight() == 0) {
			fwdResult.best_read_index = fwdResult.best_ref_index = 0;
		}

		return curr_max;
	}

	inline AlignmentMatrixFast::Score FastUnrolledfwdFillMatrixLine(char const * const refSeq,
			char const * const qrySeq, FwdResults & fwdResult, int readId,
			AlignmentMatrixFast::Score diag_score,
			AlignmentMatrixFast::MatrixElement const & up,
			AlignmentMatrixFast::MatrixElement const & left,
			AlignmentMatrixFast::MatrixElement * current,
			char & currentDirection,
			AlignmentMatrixFast::Score & curr_max,
			int x,
			int y,
			char read_char_cache) {
			
			bool const eq = read_char_cache == refSeq[x];
			AlignmentMatrixFast::Score const diag_cell = diag_score + (eq*mat) + ((!eq)*mis);

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

		return curr_max;
	}

	/**
	 * Follows direction matrix (backtracking)
	 * Creates reversed binary representation of the cigar string
	 */
	bool revBacktrack(char const * const refSeq, char const * const qrySeq,
			FwdResults & fwdResults, int readId);

	/**
	 * Converts CIGAR to text and computes MD flag + stats
	 */
	int convertCigar(char const * const refSeq, Align & result,
			FwdResults & fwdResults, int const externalQStart,
			int const externalQEnd);


	int NumberOfSetBits(uint32_t i);

	void addPosition(Align & result, int & nmIndex, int posInRef, int posInRead, int Yi);

};

#endif /* CONVEXALIGN_H_ */

}  // namespace Convex
