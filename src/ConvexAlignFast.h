/**
 * Contact: philipp.rescheneder@gmail.com
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
	ConvexAlignFast(int const stdOutMode,
			float const match,
			float const mismatch,
			float const gapOpen,
			float const gapExtend,
			float const gapExtendMin,
			float const gapDecay);

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
	 * Default Max length of binaryCigar
	 */
	int const defaultMaxBinaryCigarLength;

	/**
	 * Max length of binaryCigar
	 * If too small -> reallocate -> align -> set back to default
	 */
	int maxBinaryCigarLength;


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
	int const stdoutPrintAlignCorridor;

	/**
	 * Debug
	 */
	int alignmentId;

	/**
	 * Print CIGAR element to cigar string
	 * Helper function used by convertCigar
	 */
	int printCigarElement(char const op, int const length, char * cigar, int & cigarOpCount, int cigarMaxLength);

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

	AlignmentMatrixFast::Score fwdFillMatrixSSESimple(char const * const refSeq,
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

/*#define FILLMATRIX(DIAG_SCORE,UP_REF,LEFT_REF,CURRENT,CURRENT_DIR) {\
			diag_score=DIAG_SCORE;\
			AlignmentMatrixFast::MatrixElement const & up = UP_REF;\
			AlignmentMatrixFast::MatrixElement const & left = LEFT_REF;\
			AlignmentMatrixFast::MatrixElement * current=CURRENT;\
			char & currentDirection=CURRENT_DIR;\
			bool const eq = read_char_cache == refSeq[x]; \
			diag_cell = diag_score + (eq ? mat : mis);\
\
			up_cell = 0;\
			left_cell = 0;\
\
			ins_run = 0;\
			del_run = 0;\
\
			if (up.direction == CIGAR_I) {\
				ins_run = up.indelRun;\
				if (up.score == 0) {\
					up_cell = 0;\
				} else {\
					up_cell = up.score\
							+ std::min(gap_ext_min,\
									gap_ext + ins_run * gap_decay);\
				}\
			} else {\
				up_cell = up.score + gap_open_read;\
			}\
\
			if (left.direction == CIGAR_D) {\
				del_run = left.indelRun;\
				if (left.score == 0) {\
					left_cell = 0;\
				} else {\
					left_cell = left.score\
							+ std::min(gap_ext_min,\
									gap_ext + del_run * gap_decay);\
				}\
			} else {\
				left_cell = left.score + gap_open_ref;\
			}\
\
			AlignmentMatrixFast::Score max_cell = 0;\
			max_cell = std::max(left_cell, max_cell);\
			max_cell = std::max(diag_cell, max_cell);\
			max_cell = std::max(up_cell, max_cell);\
\
			if (del_run > 0 && max_cell == left_cell) {\
				current->score = max_cell;\
				current->direction = CIGAR_D;\
				currentDirection = CIGAR_D;\
				current->indelRun = del_run + 1;\
			} else if (ins_run > 0 && max_cell == up_cell) {\
				current->score = max_cell;\
				current->direction = CIGAR_I;\
				currentDirection = CIGAR_I;\
				current->indelRun = ins_run + 1;\
			} else if (max_cell == diag_cell) {\
				current->score = max_cell;\
				if (eq) {\
					current->direction = CIGAR_EQ;\
					currentDirection = CIGAR_EQ;\
				} else {\
					current->direction = CIGAR_X;\
					currentDirection = CIGAR_X;\
				}\
				current->indelRun = 0;\
			} else if (max_cell == left_cell) {\
				current->score = max_cell;\
				current->direction = CIGAR_D;\
				currentDirection = CIGAR_D;\
				current->indelRun = 1;\
			} else if (max_cell == up_cell) {\
				current->score = max_cell;\
				current->direction = CIGAR_I;\
				currentDirection = CIGAR_I;\
				current->indelRun = 1;\
			} else {\
				current->score = 0;\
				current->direction = CIGAR_STOP;\
				currentDirection = CIGAR_STOP;\
				current->indelRun = 0;\
			}\
\
			if (max_cell > curr_max) {\
				curr_max = max_cell;\
				fwdResult.best_ref_index = x;\
				fwdResult.best_read_index = y;\
				fwdResult.max_score = curr_max;\
			}\
		}

	inline AlignmentMatrixFast::Score FastUnrolledfwdFillMatrixMaster(char const * const refSeq,
			char const * const qrySeq, FwdResults & fwdResult, int readId) {

		AlignmentMatrixFast::Score curr_max = -1.0f;
		int const matrix_height=matrix->getHeight();


		AlignmentMatrixFast::Score diag_score;
		AlignmentMatrixFast::Score diag_cell;
		AlignmentMatrixFast::Score left_cell;
		AlignmentMatrixFast::Score up_cell;

		int ins_run;
		int del_run;


		{
			int y=0;

			matrix->prepareLine(y);
			int xOffset = matrix->getCorridorOffset(y);
			int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

			char const read_char_cache = qrySeq[y];

			for (int x = std::max(0, xOffset); x < xMax; ++x) {
				FILLMATRIX(matrix->empty.score,matrix->empty,*matrix->getElementCurr(x - 1, y),matrix->getElementEditCurrFast(x,y),*matrix->getDirectionCurrFast(x, y));
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
					FILLMATRIX(matrix->getElementUp(x - 1, y - 1)->score,*matrix->getElementUp(x,y - 1),*matrix->getElementCurr(x - 1, y),matrix->getElementEditCurrFast(x,y),*matrix->getDirectionCurrFast(x, y));
				}			
			} else {
				for(int x=xMin; x<=xFastMin; x++)
				{
					FILLMATRIX(matrix->getElementUp(x - 1, y - 1)->score,*matrix->getElementUp(x,y - 1),*matrix->getElementCurr(x - 1, y),matrix->getElementEditCurrFast(x,y),*matrix->getDirectionCurrFast(x, y));
				}

				for(int x=xFastMin+1; x<xFastMax; x++)
				{
					FILLMATRIX(matrix->getElementUpFast(x - 1, y - 1)->score,*matrix->getElementUpFast(x,y - 1),*matrix->getElementCurrFast(x - 1, y),matrix->getElementEditCurrFast(x,y),*matrix->getDirectionCurrFast(x, y));
				}

				for(int x=xFastMax; x<xMax; x++)
				{
					FILLMATRIX(matrix->getElementUp(x - 1, y - 1)->score,*matrix->getElementUp(x,y - 1),*matrix->getElementCurr(x - 1, y),matrix->getElementEditCurrFast(x,y),*matrix->getDirectionCurrFast(x, y));
				}			
			}
		}
		fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
		if (matrix->getHeight() == 0) {
			fwdResult.best_read_index = fwdResult.best_ref_index = 0;
		}

		return curr_max;
	}

	// inline AlignmentMatrixFast::Score FastUnrolledfwdFillMatrixMaster(char const * const refSeq,
	// 		char const * const qrySeq, FwdResults & fwdResult, int readId) {

	// 	AlignmentMatrixFast::Score curr_max = -1.0f;
	// 	int matrix_height=matrix->getHeight();


	// 	{
	// 		int y=0;

	// 		matrix->prepareLine(y);
	// 		int xOffset = matrix->getCorridorOffset(y);
	// 		int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

	// 		for (int x = std::max(0, xOffset); x < xMax; ++x) {
	// 			ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 				matrix->empty.score,
	// 				matrix->empty,
	// 				*matrix->getElementCurr(x - 1, y),
	// 				matrix->getElementEditCurrFast(x,y),
	// 				*matrix->getDirectionCurrFast(x, y),
	// 				curr_max,
	// 				x,
	// 				y,
	// 				qrySeq[0]);
	// 		}
	// 	}


	// 	for (int y = 1; y < matrix_height; ++y) {

	// 		matrix->prepareLine(y);
	// 		int xOffset = matrix->getCorridorOffset(y);

	// 		int xMin=std::max(0, xOffset);
	// 		int xMax=std::min(xOffset + matrix->getCorridorLength(y),matrix->getWidth());

	// 		int xFastMin=std::max(0,matrix->getCorridorOffset(y-1)+1);
	// 		int xFastMax=std::min(matrix->getCorridorOffset(y-1)+matrix->getCorridorLength(y-1),matrix->getWidth());


	// 		char const read_char_cache = qrySeq[y];

	// 		if(xMin>xFastMin || xFastMax>xMax)
	// 		{
	// 			for(int x=xMin; x<xMax; x++)
	// 			{
	// 				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 					matrix->getElementUp(x - 1, y - 1)->score,
	// 					*matrix->getElementUp(x,y - 1),
	// 					*matrix->getElementCurr(x - 1, y),
	// 					matrix->getElementEditCurrFast(x,y),
	// 					*matrix->getDirectionCurrFast(x, y),
	// 					curr_max,
	// 					x,
	// 					y,
	// 					read_char_cache);
	// 			}			
	// 		} else {
	// 			for(int x=xMin; x<=xFastMin; x++)
	// 			{
	// 				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 					matrix->getElementUp(x - 1, y - 1)->score,
	// 					*matrix->getElementUp(x,y - 1),
	// 					*matrix->getElementCurr(x - 1, y),
	// 					matrix->getElementEditCurrFast(x,y),
	// 					*matrix->getDirectionCurrFast(x, y),
	// 					curr_max,
	// 					x,
	// 					y,
	// 					read_char_cache);
	// 			}

	// 			for(int x=xFastMin+1; x<xFastMax; x++)
	// 			{
	// 				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 					matrix->getElementUpFast(x - 1, y - 1)->score,
	// 					*matrix->getElementUpFast(x,y - 1),
	// 					*matrix->getElementCurrFast(x - 1, y),
	// 					matrix->getElementEditCurrFast(x,y),
	// 					*matrix->getDirectionCurrFast(x, y),
	// 					curr_max,
	// 					x,
	// 					y,
	// 					read_char_cache);
	// 			}


	// 			for(int x=xFastMax; x<xMax; x++)
	// 			{
	// 				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 					matrix->getElementUp(x - 1, y - 1)->score,
	// 					*matrix->getElementUp(x,y - 1),
	// 					*matrix->getElementCurr(x - 1, y),
	// 					matrix->getElementEditCurrFast(x,y),
	// 					*matrix->getDirectionCurrFast(x, y),
	// 					curr_max,
	// 					x,
	// 					y,
	// 					read_char_cache);		
	// 			}			
	// 		}


	 
	// 		for (int x = std::max(0, xOffset); x < xMax; ++x) {
	// 			if(y==0 ||
	// 			  (x-1 <= (matrix->getCorridorOffset(y-1)) ) ||
	// 			  (x >= (matrix->getCorridorOffset(y-1)+matrix->getCorridorLength(y-1)) )  ) {
	// 				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 					matrix->getElementUp(x - 1, y - 1)->score,
	// 					*matrix->getElementUp(x,y - 1),
	// 					*matrix->getElementCurr(x - 1, y),
	// 					matrix->getElementEditCurrFast(x,y),
	// 					*matrix->getDirectionCurrFast(x, y),
	// 					curr_max,
	// 					x,
	// 					y);
	// 			} else {
	// 				ConvexAlignFast::FastUnrolledfwdFillMatrixLine(refSeq,qrySeq,fwdResult,readId,
	// 					matrix->getElementUpFast(x - 1, y - 1)->score,
	// 					*matrix->getElementUpFast(x,y - 1),
	// 					*matrix->getElementCurrFast(x - 1, y),
	// 					matrix->getElementEditCurrFast(x,y),
	// 					*matrix->getDirectionCurrFast(x, y),
	// 					curr_max,
	// 					x,
	// 					y);
	// 			}
	// 		}

	// 	}
	// 	fwdResult.qend = (matrix->getHeight() - fwdResult.best_read_index) - 1;
	// 	if (matrix->getHeight() == 0) {
	// 		fwdResult.best_read_index = fwdResult.best_ref_index = 0;
	// 	}

	// 	return curr_max;
	// }

	// inline void FastUnrolledfwdFillMatrixLine(char const * const refSeq,
	// 		char const * const qrySeq, FwdResults & fwdResult, int readId,
	// 		AlignmentMatrixFast::Score diag_score,
	// 		AlignmentMatrixFast::MatrixElement const & up,
	// 		AlignmentMatrixFast::MatrixElement const & left,
	// 		AlignmentMatrixFast::MatrixElement * current,
	// 		char & currentDirection,
	// 		AlignmentMatrixFast::Score & curr_max,
	// 		int x,
	// 		int y,
	// 		char read_char_cache) {
			
	// 		bool const eq = read_char_cache == refSeq[x];
	// 		AlignmentMatrixFast::Score const diag_cell = diag_score + (eq ? mat : mis);

	// 		AlignmentMatrixFast::Score up_cell = 0;
	// 		AlignmentMatrixFast::Score left_cell = 0;

	// 		int ins_run = 0;
	// 		int del_run = 0;

	// 		if (up.direction == CIGAR_I) {
	// 			ins_run = up.indelRun;
	// 			if (up.score == 0) {
	// 				up_cell = 0;
	// 			} else {
	// 				up_cell = up.score
	// 						+ std::min(gap_ext_min,
	// 								gap_ext + ins_run * gap_decay);
	// 			}
	// 		} else {
	// 			up_cell = up.score + gap_open_read;
	// 		}

	// 		if (left.direction == CIGAR_D) {
	// 			del_run = left.indelRun;
	// 			if (left.score == 0) {
	// 				left_cell = 0;
	// 			} else {
	// 				left_cell = left.score
	// 						+ std::min(gap_ext_min,
	// 								gap_ext + del_run * gap_decay);
	// 			}
	// 		} else {
	// 			left_cell = left.score + gap_open_ref;
	// 		}

	// 		//find max
	// 		AlignmentMatrixFast::Score max_cell = 0;
	// 		max_cell = std::max(left_cell, max_cell);
	// 		max_cell = std::max(diag_cell, max_cell);
	// 		max_cell = std::max(up_cell, max_cell);

	// 		int a=del_run > 0 && max_cell == left_cell;
	// 		int b=(!a) && ins_run > 0 && max_cell == up_cell;
	// 		int c=(!a && !b) && max_cell == diag_cell;
	// 		int d=(!a && !b && !c) && max_cell == left_cell;
	// 		int e=(!a && !b && !c && !d) && max_cell == up_cell;
	// 		int f=!(a*b*c*d*e);

	// 		current->score = a*max_cell;
	// 		current->direction = a*CIGAR_D;
	// 		currentDirection = a*CIGAR_D;
	// 		current->indelRun = a*(del_run + 1);

	// 		current->score = (!b)*current->score + b*max_cell;
	// 		current->direction = (!b)*current->direction + b*CIGAR_I;
	// 		currentDirection = (!b)*currentDirection + b*CIGAR_I;
	// 		current->indelRun = (!b)*current->indelRun + b*(ins_run + 1);

	// 		current->score = (!c)*current->score + c*max_cell;
	// 		current->direction = (!c)*current->direction + (c*eq*CIGAR_EQ) + (c*(!eq)*CIGAR_X);
	// 		currentDirection = (!c)*currentDirection + (c*eq*CIGAR_EQ) + (c*(!eq)*CIGAR_X);
	// 		current->indelRun = (!c)*current->indelRun + c*0;

	// 		current->score = (!d)*current->score + d*max_cell;
	// 		current->direction = (!d)*current->direction + d*CIGAR_D;
	// 		currentDirection = (!d)*currentDirection + d*CIGAR_D;
	// 		current->indelRun = (!d)*current->indelRun + d*1;

	// 		current->score = (!e)*current->score + e*max_cell;
	// 		current->direction = (!e)*current->direction + e*CIGAR_I;
	// 		currentDirection = (!e)*currentDirection + e*CIGAR_I;
	// 		current->indelRun = (!e)*current->indelRun + e*1;

	// 		current->score = (!f)*current->score + 0;
	// 		current->direction = (!f)*current->direction + f*CIGAR_STOP;
	// 		currentDirection = (!f)*currentDirection + f*CIGAR_STOP;
	// 		current->indelRun = (!f)*current->indelRun + 0;*/

	// 		if (del_run > 0 && max_cell == left_cell) {
	// 			current->score = max_cell;
	// 			current->direction = CIGAR_D;
	// 			currentDirection = CIGAR_D;
	// 			current->indelRun = del_run + 1;
	// 		} else if (ins_run > 0 && max_cell == up_cell) {
	// 			current->score = max_cell;
	// 			current->direction = CIGAR_I;
	// 			currentDirection = CIGAR_I;
	// 			current->indelRun = ins_run + 1;
	// 		} else if (max_cell == diag_cell) {
	// 			current->score = max_cell;
	// 			if (eq) {
	// 				current->direction = CIGAR_EQ;
	// 				currentDirection = CIGAR_EQ;
	// 			} else {
	// 				current->direction = CIGAR_X;
	// 				currentDirection = CIGAR_X;
	// 			}
	// 			current->indelRun = 0;
	// 		} else if (max_cell == left_cell) {
	// 			current->score = max_cell;
	// 			current->direction = CIGAR_D;
	// 			currentDirection = CIGAR_D;
	// 			current->indelRun = 1;
	// 		} else if (max_cell == up_cell) {
	// 			current->score = max_cell;
	// 			current->direction = CIGAR_I;
	// 			currentDirection = CIGAR_I;
	// 			current->indelRun = 1;
	// 		} else {
	// 			current->score = 0;
	// 			current->direction = CIGAR_STOP;
	// 			currentDirection = CIGAR_STOP;
	// 			current->indelRun = 0;
	// 		}

	// 		if (max_cell > curr_max) {
	// 			curr_max = max_cell;
	// 			fwdResult.best_ref_index = x;
	// 			fwdResult.best_read_index = y;
	// 			fwdResult.max_score = curr_max;
	// 		}
	// }*/

	/**
	 * Follows direction matrix (backtracking)
	 * Creates reversed binary representation of the cigar string
	 */
	bool revBacktrack(char const * const refSeq, char const * const qrySeq,
			FwdResults & fwdResults, int readId);

	/**
	 * Converts CIGAR to text and computes MD flag + stats
	 */
	int convertCigar(char const * const refSeq, int const refSeqLength, Align & result,
			FwdResults & fwdResults, int const externalQStart,
			int const externalQEnd);


	int NumberOfSetBits(uint32_t i);

	void addPosition(Align & result, int & nmIndex, int posInRef, int posInRead, int Yi);
	void checkMdBufferLength(int md_offset, Align& result, int const minDiff = 0);
};

#endif /* CONVEXALIGN_H_ */

}  // namespace Convex
