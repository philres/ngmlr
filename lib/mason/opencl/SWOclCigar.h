/*
 * SWOclCigar.h
 *
 *  Created on: May 25, 2011
 *      Author: philipp_
 */

#ifndef SWOCLCIGAR_H_
#define SWOCLCIGAR_H_

#include "SWOcl.h"

#undef module_name
#define module_name "Cigar (OpenCL)"

class SWOclCigar: public SWOcl {
public:
	SWOclCigar(OclHost * host);
	virtual ~SWOclCigar();
	virtual int BatchAlign(int const mode, int const batchSize, char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results, void * extData);

	virtual int GetAlignBatchSize() const;

private:

	cl_kernel swAlignScoreKernel;
	cl_kernel swAlignScoreKernelGlobal;
	cl_kernel swAlignBacktrackingKernel;

	int batch_size_align;

	void runSwBatchKernel(cl_kernel swScoreAlign, const int batchSize, const char * const * const qrySeqList,
			const char * const * const refSeqList, char * bsDirection, cl_mem & results_gpu, cl_mem & alignments, short * const result,
			short * calignments, cl_mem & matrix_gpu, cl_mem & bsdirection_gpu);
	bool computeCigarMD(Align & result, int const gpuCigarOffset, short const * const gpuCigar, char const * const refSeq, char const * const qrySeq, char const bsFrom, char const bsTo);

	int computeAlignmentBatchSize();

	long getMaxAllocSize(int const batch_size);
};

#endif /* SWOCLCIGAR_H_ */
