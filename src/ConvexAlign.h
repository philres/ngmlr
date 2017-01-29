/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef CONVEXALIGN_H_
#define CONVEXALIGN_H_

#include "IAlignment.h"

#include "Types.h"
#include "AlignmentMatrix.h"

namespace Convex {

class ConvexAlign: public IAlignment {

public:
	ConvexAlign(int const stdOutMode,
			float const match,
			float const mismatch,
			float const gapOpen,
			float const gapExtend,
			float const gapExtendMin,
			float const gapDecay);

	virtual ~ConvexAlign();

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
	AlignmentMatrix::Score mat;
	AlignmentMatrix::Score mis;
	AlignmentMatrix::Score gap_open_read;
	AlignmentMatrix::Score gap_open_ref;
	AlignmentMatrix::Score gap_ext;
	AlignmentMatrix::Score gap_ext_min;
	AlignmentMatrix::Score gap_decay;

	AlignmentMatrix * matrix;

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
	int const stdoutPrintAlignCorridor;

	/**
	 * Print CIGAR element to cigar string
	 * Helper function used by convertCigar
	 */
	int printCigarElement(char const op, int const length, char * cigar);

	/**
	 * Forward step of alignment
	 * Fills directionMatrix
	 */
	AlignmentMatrix::Score fwdFillMatrix(char const * const refSeq,
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



};

#endif /* CONVEXALIGN_H_ */

}  // namespace Convex
