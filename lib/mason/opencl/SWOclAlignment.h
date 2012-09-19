/*
 * SWOclAlignment.h
 *
 *  Created on: May 26, 2011
 *      Author: philipp_
 */

#ifndef SWOCLALIGNMENT_H_
#define SWOCLALIGNMENT_H_

#include "SWOcl.h"

#undef module_name
#define module_name "Alignment (OpenCL)"

class SWOclAlignment : public SWOcl {
public:
	SWOclAlignment(OclHost * host);
	virtual ~SWOclAlignment();

	//virtual int BatchAlign(const int mode, const int batchSize, const char * const * const refSeqList, const char * const * const qrySeqList, Align * const results,void * extData);
	virtual int BatchAlign(const int mode, const int batchSize, const char * const * const refSeqList, const char * const * const qrySeqList, Align * const results,void * extData);

	virtual int GetAlignBatchSize() const;
private:

	cl_kernel swAlignScoreKernel;
	cl_kernel swAlignBacktrackingKernel;
	cl_kernel swAlignScoreKernelGlobal;

	int batch_size_align;

	void runSwBatchKernel(cl_kernel swScoreAlign, const int batchSize, const char * const * const qrySeqList, const char * const * const refSeqList, cl_mem & results_gpu, cl_mem & alignments, short * const result, char * calignments, cl_mem & matrix_gpu);
	int computeAlignmentBatchSize();

	long getMaxAllocSize(int const batch_size);
};

#endif /* SWOCLALIGNMENT_H_ */
