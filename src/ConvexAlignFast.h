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
