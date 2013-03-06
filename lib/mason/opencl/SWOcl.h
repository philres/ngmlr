/*
 * SWOcl.h
 *
 *  Created on: Apr 8, 2011
 *      Author: philipp_
 */

#ifndef SWOCL_H_
#define SWOCL_H_

#include <CL/opencl.h>

#include "IAlignment.h"
#include "ILog.h"
#include "IConfig.h"

//#define _DEBUG
//#define _TEST

#undef module_name
#define module_name "Score (OpenCL)"

class OclHost;

#define CIGAR_M 0
#define CIGAR_I 1
#define CIGAR_D 2
#define CIGAR_N 3
#define CIGAR_S 4
#define CIGAR_H 5
#define CIGAR_P 6
#define CIGAR_EQ 7
#define CIGAR_X 8

class SWOcl: public IAlignment {
public:
	//SWOcl(const int gpu_id);
	SWOcl(char const * const oclSwScoreSourceCode, char const * const additional_defines, OclHost * host);
	virtual ~SWOcl();
	virtual int GetScoreBatchSize() const;
	virtual int GetAlignBatchSize() const = 0;
	virtual int BatchScore(const int mode, const int batchSize, const char * const * const refSeqList,
			const char * const * const qrySeqList, float * const results, void * extData);
	virtual int BatchAlign(const int mode, const int batchSize, const char * const * const refSeqList,
			const char * const * const qrySeqList, Align * const results, void * extData) = 0;

	OclHost * getHost() {
		return host;
	}

private:
	//unsigned int block_count;
	int batch_size_scoring;
	//int batch_size_align;

protected:

	/*Number of threads executed by each multi processor*/
	static unsigned int const threads_per_block = 256;
	unsigned int block_multiplier;
	unsigned int step_count;

	int alignment_length;
	int matrix_length;
	int config_ref_size;

	OclHost * host;
	//static int programUserCount;
//	static cl_program clProgram;
	cl_program clProgram;

	cl_kernel swScoringKernel;
	cl_kernel swScoringKernelGlobal;
	cl_kernel interleaveKernel;

	cl_mem c_scaff_gpu;
	cl_mem c_reads_gpu;
	cl_mem reads_gpu;
	cl_mem scaffold_gpu;
	cl_mem cmPinnedBufRead;
	cl_mem cmPinnedBufRef;

	int read_data_size;
	char * cpu_read_data;
	int ref_data_size;
	char * cpu_ref_data;

	int alignments_per_thread;

	void checkMemory();

	void copySeqDataToDevice(char *cpu_read_data, char *cpu_ref_data, const char * const * const reads, const char * const * const refs,
			const int alignment_number, const int batch_size);
	void copy(char const * const * const refs, char * cpu_ref_data, int const alignment_number, int const batch_size, int const refsize);
	int runSwScoreKernel(cl_kernel scoreKernel, const int batchSize, const char * const * const qrySeqList,
			const char * const * const refSeqList, char * bsDirection, cl_mem & bsdirection_gpu, cl_mem & results_gpu,
			float * const results);
	//void copySeqDataToDeviceAlign(char *cpu_read_data, char *cpu_ref_data, const char * const * const reads, const char * const * const refs, const int alignment_number, const int batch_size);

	virtual int computeAlignmentBatchSize() = 0;

private:
	unsigned int computeScoringBatchSize();
	int getMaxAllocSize();

};

#endif /* SWOCL_H_ */
