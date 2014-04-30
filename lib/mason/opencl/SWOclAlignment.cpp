/*
 * SWOclAlignment.cpp
 *
 *  Created on: May 26, 2011
 *      Author: philipp_
 */

#include "SWOclAlignment.h"

#include <stdio.h>
#include <string.h>
#include <sstream>

#include "Timing.h"
#include "OclHost.h"

#include "oclSwAlignment.h"

#define pRef pBuffer1
#define pQry pBuffer2

int const result_number = 4;

SWOclAlignment::SWOclAlignment(OclHost * host) :
		SWOcl(
				oclSwAlignment,
				" -D result_number=4 -D CIGAR_M=0 -D CIGAR_I=1 -D CIGAR_D=2 -D CIGAR_N=3 -D CIGAR_S=4 -D CIGAR_H=5 -D CIGAR_P=6 -D CIGAR_EQ=7 -D CIGAR_X=8 ",
				host) {
	batch_size_align = computeAlignmentBatchSize();
	swAlignScoreKernel = host->setupKernel(clProgram, "oclSW_Score");
	swAlignScoreKernelGlobal = host->setupKernel(clProgram,
			"oclSW_ScoreGlobal");
	swAlignBacktrackingKernel = host->setupKernel(clProgram,
			"oclSW_Backtracking");
}

SWOclAlignment::~SWOclAlignment() {
	clReleaseKernel(swAlignBacktrackingKernel);
	clReleaseKernel(swAlignScoreKernel);
	clReleaseKernel(swAlignScoreKernelGlobal);
}

//#define _TEST

#ifdef _TEST
#include "Timer.h"

#define benchCPP(x,z)  test.start(); for(int y=0;y<100;++y) x; test.stop(); Log.Message("%s\t%f", z, test.getElapsedTimeInMilliSec());
#define benchOCL(x,z)  test.start(); for(int y=0;y<100;++y) x;host->waitForDevice(); test.stop(); Log.Message("%s\t%f", z, test.getElapsedTimeInMilliSec());
#else
#define benchCPP(x,z) x;
#define benchOCL(x,z) x;
#endif

void SWOclAlignment::runSwBatchKernel(cl_kernel swScoreAlign,
		const int batchSize, const char * const * const qrySeqList,
		const char * const * const refSeqList, cl_mem & results_gpu,
		cl_mem & alignments, short * const result, char * calignments,
		cl_mem & matrix_gpu) {
#ifdef _TEST
	Timer2 test;
#endif

	const size_t cnDim = batch_size_align / alignments_per_thread;
	const size_t cBlockSize = threads_per_block;
	int runbatchSize = std::min(batch_size_align, batchSize);
	copySeqDataToDevice(cpu_read_data, cpu_ref_data, qrySeqList, refSeqList,
			runbatchSize, batch_size_align);
	for (int i = 0; i < batchSize; i += batch_size_align) {

		if (host->isGPU()) {
			host->writeToDevice(scaffold_gpu, CL_FALSE, 0,
					ref_data_size * sizeof(cl_char), cpu_ref_data);
			host->writeToDevice(reads_gpu, CL_FALSE, 0,
					read_data_size * sizeof(cl_char), cpu_read_data);
			host->waitForDevice();
			benchOCL(host->executeKernel(interleaveKernel, cnDim, cBlockSize),
					"Intel:");
		}

		benchOCL(host->executeKernel(swScoreAlign, cnDim, cBlockSize),
				"Score:");

		benchOCL(
				host->executeKernel(swAlignBacktrackingKernel, cnDim, cBlockSize),
				"Backtr:");

		int nextBatch = (i + batch_size_align);
		int nextRunbatchSize = std::min(batch_size_align,
				batchSize - nextBatch);
		if (nextRunbatchSize > 0) {
			copySeqDataToDevice(cpu_read_data, cpu_ref_data,
					qrySeqList + nextBatch, refSeqList + nextBatch,
					nextRunbatchSize, batch_size_align);
		}
		host->readFromDevice(results_gpu, CL_FALSE, 0,
				result_number * runbatchSize, result + i * result_number,
				sizeof(cl_short));
		host->readFromDevice(alignments, CL_FALSE, 0,
				runbatchSize * alignment_length * 2,
				calignments + alignment_length * 2 * i, sizeof(cl_char));
		runbatchSize = nextRunbatchSize;

	}
	host->waitForDevice();
}

int SWOclAlignment::BatchAlign(int const mode, int const batchSize_,
		char const * const * const refSeqList_,
		char const * const * const qrySeqList_, Align * const results,
		void * extData) {

	if (batchSize_ <= 0) {
		Log.Warning("Align for batchSize <= 0");
		return 0;
	}

	bool batchSizeDif = !host->isGPU() && (batchSize_ % 4 != 0);
//	int batchSize = (batchSizeDif) ? batchSize_ + 4 : batchSize_;
	int batchSize = (batchSizeDif) ? ((int)(batchSize_ / 4) + 1) * 4 : batchSize_;
	char const * * const tmpRefSeqList = new char const *[batchSize];
	char const * * const tmpQrySeqList = new char const *[batchSize];

	if (batchSizeDif) {
		Log.Warning("Batchsize not a multiple of 4.");
		for (int i = 0; i < batchSize_; ++i) {
			tmpRefSeqList[i] = refSeqList_[i];
			tmpQrySeqList[i] = qrySeqList_[i];
		}
		for (int i = batchSize_; i < batchSize; ++i) {
			tmpRefSeqList[i] = refSeqList_[0];
			tmpQrySeqList[i] = qrySeqList_[0];
		}
	}
	char const * const * const refSeqList =
			batchSizeDif ? tmpRefSeqList : refSeqList_;
	char const * const * const qrySeqList =
			batchSizeDif ? tmpQrySeqList : qrySeqList_;

	cl_kernel scoreKernel;
	switch ((mode & 0xFF)) {
	case 0:
		Log.Verbose("Alignment mode: local");
		scoreKernel = swAlignScoreKernel;
		break;
	case 1:
		Log.Verbose("Alignment mode: end-free");
		scoreKernel = swAlignScoreKernelGlobal;
		break;
	default:
		Log.Error("Unsupported alignment mode %i", mode & 0xFF);
		exit(-1);
	}

	Timer timer;
	timer.ST();

	cl_int ciErrNum = 0;
	cl_mem c_scaff_gpu = scaffold_gpu;
	if (host->isGPU()) {
		c_scaff_gpu = host->allocate(CL_MEM_READ_WRITE,
				ref_data_size * sizeof(cl_char));
		ciErrNum |= clSetKernelArg(interleaveKernel, 0, sizeof(cl_mem),
				(void *) (&scaffold_gpu));
		ciErrNum |= clSetKernelArg(interleaveKernel, 1, sizeof(cl_mem),
				&c_scaff_gpu);
	} else {
		c_scaff_gpu = scaffold_gpu;
	}
	cl_mem results_gpu = host->allocate(CL_MEM_READ_WRITE,
			result_number * batch_size_align * sizeof(cl_short));
	cl_mem matrix_gpu = host->allocate(
			CL_MEM_READ_WRITE,
			batch_size_align * (Config.GetInt("corridor") + 2)
					* (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char));
	cl_mem alignments_gpu = host->allocate(CL_MEM_READ_WRITE,
			batch_size_align * alignment_length * 2 * sizeof(cl_char));

	//Set parameter
	ciErrNum |= clSetKernelArg(scoreKernel, 0, sizeof(cl_mem),
			(void *) (&c_scaff_gpu));
	ciErrNum |= clSetKernelArg(scoreKernel, 1, sizeof(cl_mem),
			(void*) (&reads_gpu));
	ciErrNum |= clSetKernelArg(scoreKernel, 2, sizeof(cl_mem), &results_gpu);
	ciErrNum |= clSetKernelArg(scoreKernel, 3, sizeof(cl_mem), &matrix_gpu);

	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 0, sizeof(cl_mem),
			(void *) (&c_scaff_gpu));
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 1, sizeof(cl_mem),
			(void*) (&reads_gpu));
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 2, sizeof(cl_mem),
			&results_gpu);
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 3, sizeof(cl_mem),
			&matrix_gpu);
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 4, sizeof(cl_mem),
			&alignments_gpu);

	host->checkClError("Unable to set kernel parameters", ciErrNum);

	char * calignments = new char[batchSize * alignment_length * 2];
	short * ref_offset = new short[result_number * batchSize];

	runSwBatchKernel(scoreKernel, batchSize, qrySeqList, refSeqList,
			results_gpu, alignments_gpu, ref_offset, calignments, matrix_gpu);

	for (int i = 0; i < batchSize_; ++i) {
		//results[i].pQry = new char[alignment_length];
		//results[i].pRef = new char[alignment_length];
		float total = 0.0f;
		float match = 0.0f;
		char * read = results[i].pQry;
		char * ref = results[i].pRef;
		char * tempAlign = calignments + i * alignment_length * 2;

		if (host->isGPU()) {
			int index = 0;
			for (int j = ref_offset[result_number * i + 3] + 1;
					j < alignment_length; ++j) {
				ref[index] = tempAlign[j];
				read[index] = tempAlign[j + alignment_length];

				total += 1;
				if (read[index] == ref[index]) {
					match += 1;
				}
				index += 1;
			}
			results[i].PositionOffset = ref_offset[result_number * i];
			results[i].QStart = ref_offset[result_number * i + 1];
			results[i].QEnd = ref_offset[result_number * i + 2];
			results[i].Identity = match * 1.0f / total;
		} else {
			int index = 0;
			int offset = (i / alignments_per_thread) * result_number
					* alignments_per_thread;
			int k = (i) % alignments_per_thread;
			for (int j = ref_offset[offset + alignments_per_thread * 3 + k] + 1;
					j < alignment_length; ++j) {
				ref[index] = tempAlign[j];
				read[index] = tempAlign[j + alignment_length];

				total += 1;
				if (read[index] == ref[index]) {
					match += 1;
				}
				index += 1;
			}
			results[i].PositionOffset = ref_offset[offset + k];
			results[i].QStart = ref_offset[offset + alignments_per_thread * 1
					+ k];
			results[i].QEnd =
					ref_offset[offset + alignments_per_thread * 2 + k];
			results[i].Identity = match * 1.0f / total;
		}
	}

	delete[] calignments;
	calignments = 0;
	delete[] ref_offset;
	ref_offset = 0;

	host->waitForDevice();

	cl_int clErr = clReleaseMemObject(matrix_gpu);
	clErr = clReleaseMemObject(results_gpu);
	clErr = clReleaseMemObject(alignments_gpu);
	if (host->isGPU()) {
		clErr = clReleaseMemObject(c_scaff_gpu);
	}
	host->checkClError("Unable to release memory for alignment.", clErr);

	Log.Verbose("SW finished computing alignment for %d sequences (elapsed: %.3fs)", batchSize_, timer.ET());

	delete[] tmpRefSeqList;
	delete[] tmpQrySeqList;

	return batchSize_;

}

long SWOclAlignment::getMaxAllocSize(int const batch_size) {
	//	Log.Message("Batch size:\t %d", batch_size);
	//	Log.Message("Results:\t %d", result_number * batch_size * sizeof(cl_short));
	//
	//	Log.Message("%d %d %d %d", batch_size, (Config.GetInt("corridor") + 2), (Config.GetInt("qry_max_len") + 1), sizeof(cl_char));
	//	Log.Message("Matrix:\t %d", );
	//	Log.Message("Alignments:\t %d\n\n", );

	long corridor = (Config.GetInt("corridor") + 2);
	long qry_len = (Config.GetInt("qry_max_len") + 1);
	long s_char = sizeof(cl_char);

	long matrix = (long) batch_size * corridor * qry_len * s_char;
	//	std::cout << corridor << std::endl;
	//	std::cout << qry_len << std::endl;
	//	std::cout << s_char << std::endl;
	//	std::cout << batch_size << std::endl;
	//
	//	std::cout << std::endl << (corridor * s_char) << std::endl << (corridor * s_char * qry_len) <<
	//			std::endl << (corridor * s_char * qry_len * 92160) << std::endl;
	//	std::cout << "Matrix: " << matrix << std::endl;

	//std::cout << batch_size * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char) << std::endl;

	long alignments = (long) batch_size * (long) alignment_length * (long) 2
			* (long) sizeof(cl_char);
	//	std::cout << "alignments: " << alignments << std::endl;
	return std::max(matrix, alignments);
}

int SWOclAlignment::GetAlignBatchSize() const {
	return batch_size_align * step_count;
}

int SWOclAlignment::computeAlignmentBatchSize() {

	cl_uint mpCount = host->getDeviceInfoInt(CL_DEVICE_MAX_COMPUTE_UNITS);
	//TODO: Fix
	//	unsigned long max_alloc = host->getDeviceInfoLong(CL_DEVICE_MAX_MEM_ALLOC_SIZE) * 2.0;
	int block_count = mpCount * block_multiplier
			* (host->getThreadPerMulti() / threads_per_block);
	block_count = (block_count / mpCount) * mpCount;

	unsigned long largest_alloc = getMaxAllocSize(
			block_count * threads_per_block);
	cl_mem testAlloc = 0;
	while (!host->testAllocate(largest_alloc)) {
		block_count -= mpCount;
		largest_alloc = getMaxAllocSize(block_count * threads_per_block);
		Log.Verbose("Reducing batch size to %d", block_count * threads_per_block);
	}
//	if (!host->isGPU()) {
//		block_count *= 4;
//	}
	Log.Verbose("Multi processor count: %d", mpCount);
	Log.Verbose("Max. threads per multi processor: %d", host->getThreadPerMulti());
	Log.Verbose("Threads per block used: %d", threads_per_block);
	Log.Verbose("Block number: %d", block_count);
	Log.Verbose("Batch size: %d", (block_count * threads_per_block));
	//TODO: Print debug info

	return block_count * threads_per_block;

}
