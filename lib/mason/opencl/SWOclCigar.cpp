/*
 * SWOclCigar.cpp
 *
 *  Created on: May 25, 2011
 *      Author: philipp_
 */

#include "SWOclCigar.h"

#include <stdio.h>
#include <string.h>
#include <sstream>

#include "Timing.h"
#include "OclHost.h"

using std::stringstream;

#include "oclSwCigar.h"

#define pCigar pBuffer1
#define pMD pBuffer2

int const result_number = 4;

SWOclCigar::SWOclCigar(OclHost * host) :
		SWOcl(oclSwCigar,
				" -D result_number=4 -D CIGAR_M=0 -D CIGAR_I=1 -D CIGAR_D=2 -D CIGAR_N=3 -D CIGAR_S=4 -D CIGAR_H=5 -D CIGAR_P=6 -D CIGAR_EQ=7 -D CIGAR_X=8 ",
				host) {
	batch_size_align = computeAlignmentBatchSize();
	swAlignScoreKernel = host->setupKernel(clProgram, "oclSW_Score");
	swAlignScoreKernelGlobal = host->setupKernel(clProgram, "oclSW_ScoreGlobal");
	swAlignBacktrackingKernel = host->setupKernel(clProgram, "oclSW_Backtracking");
}

SWOclCigar::~SWOclCigar() {
	clReleaseKernel(swAlignBacktrackingKernel);
	clReleaseKernel(swAlignScoreKernel);
	clReleaseKernel(swAlignScoreKernelGlobal);
}

void SWOclCigar::runSwBatchKernel(cl_kernel swScoreAlign, const int batchSize, const char * const * const qrySeqList,
		const char * const * const refSeqList, char * bsDirection, cl_mem & results_gpu, cl_mem & alignments, short * const result,
		short * calignments, cl_mem & matrix_gpu, cl_mem & bsdirection_gpu) {
	const size_t cnDim = batch_size_align / alignments_per_thread;
	const size_t cBlockSize = threads_per_block;
	int runbatchSize = std::min(batch_size_align, batchSize);
	copySeqDataToDevice(cpu_read_data, cpu_ref_data, qrySeqList, refSeqList, runbatchSize, batch_size_align);
	for (int i = 0; i < batchSize; i += batch_size_align) {

		if (host->isGPU()) {
			host->writeToDevice(scaffold_gpu, CL_FALSE, 0, ref_data_size * sizeof(cl_char), cpu_ref_data);
			host->writeToDevice(reads_gpu, CL_FALSE, 0, read_data_size * sizeof(cl_char), cpu_read_data);

			if (bsdirection_gpu != 0) {
				host->writeToDevice(bsdirection_gpu, CL_FALSE, 0, std::min(runbatchSize, batchSize - i) * sizeof(cl_char), bsDirection + i);
			}

			host->waitForDevice();
			host->executeKernel(interleaveKernel, cnDim, cBlockSize);
		}

		host->executeKernel(swScoreAlign, cnDim, cBlockSize);
		host->executeKernel(swAlignBacktrackingKernel, cnDim, cBlockSize);

		int nextBatch = (i + batch_size_align);
		int nextRunbatchSize = std::min(batch_size_align, batchSize - nextBatch);
		if (nextRunbatchSize > 0) {
			copySeqDataToDevice(cpu_read_data, cpu_ref_data, qrySeqList + nextBatch, refSeqList + nextBatch, nextRunbatchSize,
					batch_size_align);
		}
		host->readFromDevice(results_gpu, CL_FALSE, 0, result_number * runbatchSize, result + i * result_number, sizeof(cl_short));
		host->readFromDevice(alignments, CL_FALSE, 0, runbatchSize * alignment_length * 2, calignments + alignment_length * 2 * i,
				sizeof(cl_short));
		runbatchSize = nextRunbatchSize;

	}
	host->waitForDevice();
}

int SWOclCigar::BatchAlign(int const mode, int const batchSize_, char const * const * const refSeqList_,
		char const * const * const qrySeqList_, Align * const results, void * extData) {
	if (batchSize_ <= 0) {
		Log.Warning("Align for batchSize <= 0");
		return 0;
	}

	bool batchSizeDif = !host->isGPU() && (batchSize_ % 4 != 0);
	int batchSize = (batchSizeDif) ? batchSize_ + 4 : batchSize_;
	char const * * const tmpRefSeqList = new char const *[batchSize];
	char const * * const tmpQrySeqList = new char const *[batchSize];

	if (batchSizeDif) {
		for (int i = 0; i < batchSize_; ++i) {
			tmpRefSeqList[i] = refSeqList_[i];
			tmpQrySeqList[i] = qrySeqList_[i];
		}
		for (int i = batchSize_; i < batchSize; ++i) {
			tmpRefSeqList[i] = refSeqList_[0];
			tmpQrySeqList[i] = qrySeqList_[0];
		}
	}
	char const * const * const refSeqList = batchSizeDif ? tmpRefSeqList : refSeqList_;
	char const * const * const qrySeqList = batchSizeDif ? tmpQrySeqList : qrySeqList_;

//	//	char * const * const qrySeqList = (char * const * const)qrySeqList_;
//	char * * const qrySeqList = new char *[batchSize * 2];
//
//	long corridor = (Config.GetInt("corridor"));
//	long qry_len = (Config.GetInt("qry_max_len"));
//
//	for (int i = 0; i < batchSize; ++i) {
//		std::cout << "Read: ";
//		for (int j = 0; j < (qry_len); ++j) {
//			if (qrySeqList[i][j] != '\0') {
//				std::cout << qrySeqList[i][j];
//			} else {
//				std::cout << "X";
//			}
//		}
//		std::cout << std::endl;
//		std::cout << "Ref:  ";
//		for (int j = 0; j < (corridor + qry_len); ++j) {
//			if (refSeqList[i][j] != '\0') {
//				std::cout << refSeqList[i][j];
//			} else {
//				std::cout << "X";
//			}
//		}
//		std::cout << std::endl;
//		//		std::cout << std::endl;
//		//		for (int j = 0; j < (corridor + qry_len); ++j) {
//		//			std::cout << refSeqList[i][j];
//		//		}
//		//		std::cout << std::endl;
//	}
//	std::cout << "===========================================================" << std::endl;
//
//	bool found = false;
//	for (int i = 0; i < batchSize; ++i) {
//		qrySeqList[i] = new char[qry_len + 1];
//		found = false;
//		for (int j = 0; j < (qry_len); ++j) {
//			found = found || (qrySeqList_[i][j] == '\0');
//			if (found) {
//				qrySeqList[i][j] = '\0';
//			} else {
//				qrySeqList[i][j] = qrySeqList_[i][j];
//			}
//			//			std::cout << qrySeqList[i][j];
//		}
//		//		std::cout << std::endl;
//		//		for (int j = 0; j < (corridor + qry_len); ++j) {
//		//			std::cout << refSeqList[i][j];
//		//		}
//		//		std::cout << std::endl;
//	}

//	Log.Error("Batch align: %d %d", batchSize, Config.GetInt("qry_max_len"));

	cl_kernel scoreKernel;
	switch ((mode & 0xFF)) {
	case 0:
//#ifndef NDEBUG
		Log.Verbose("Alignment mode: local");
//#endif
		scoreKernel = swAlignScoreKernel;
		break;
	case 1:
//#ifndef NDEBUG
		Log.Verbose("Alignment mode: end-free");
//#endif
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
		c_scaff_gpu = host->allocate(CL_MEM_READ_WRITE, ref_data_size * sizeof(cl_char));
		ciErrNum |= clSetKernelArg(interleaveKernel, 0, sizeof(cl_mem), (void *) (&scaffold_gpu));
		ciErrNum |= clSetKernelArg(interleaveKernel, 1, sizeof(cl_mem), &c_scaff_gpu);
	} else {
		c_scaff_gpu = scaffold_gpu;
	}

	cl_mem results_gpu = host->allocate(CL_MEM_READ_WRITE, result_number * batch_size_align * sizeof(cl_short));
	cl_mem matrix_gpu = host->allocate(CL_MEM_READ_WRITE,
			batch_size_align * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char));

	cl_mem alignments_gpu = host->allocate(CL_MEM_READ_WRITE, batch_size_align * alignment_length * 2 * sizeof(cl_short));

	//Set parameter
	ciErrNum |= clSetKernelArg(scoreKernel, 0, sizeof(cl_mem), (void *) (&c_scaff_gpu));
	ciErrNum |= clSetKernelArg(scoreKernel, 1, sizeof(cl_mem), (void *) (&reads_gpu));
	ciErrNum |= clSetKernelArg(scoreKernel, 2, sizeof(cl_mem), &results_gpu);
	ciErrNum |= clSetKernelArg(scoreKernel, 3, sizeof(cl_mem), &matrix_gpu);
	cl_mem bsdirection_gpu = 0;
	static bool const bsMapping = Config.GetInt("bs_mapping") == 1;
	if (bsMapping) {
		bsdirection_gpu = host->allocate(CL_MEM_READ_ONLY, batch_size_align * sizeof(cl_char));
		ciErrNum |= clSetKernelArg(scoreKernel, 4, sizeof(cl_mem), (void *) (&bsdirection_gpu));
	}

	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 0, sizeof(cl_mem), (void *) (&c_scaff_gpu));
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 1, sizeof(cl_mem), (void*) (&reads_gpu));
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 2, sizeof(cl_mem), &results_gpu);
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 3, sizeof(cl_mem), &matrix_gpu);
	ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 4, sizeof(cl_mem), &alignments_gpu);

	host->checkClError("Unable to set kernel parameters", ciErrNum);

	short * calignments = new short[batchSize * alignment_length * 2];
	short * gpu_return_values = new short[result_number * batchSize];

	Timer gpuTimer;
	gpuTimer.ST();
	runSwBatchKernel(scoreKernel, batchSize, qrySeqList, refSeqList, (char *) extData, results_gpu, alignments_gpu, gpu_return_values,
			calignments, matrix_gpu, bsdirection_gpu);

	Log.Verbose("GPU Time: %.3fs", gpuTimer.ET());

	Timer cpuTimer;
	cpuTimer.ST();
	for (int i = 0; i < batchSize_; ++i) {
		short * gpuCigar = calignments + i * alignment_length * 2;
		//results[i].pCigar = new char[alignment_length];
		//results[i].pMD = new char[alignment_length];
		int offset = (i / alignments_per_thread) * result_number * alignments_per_thread;
		int k = (i) % alignments_per_thread;

		char bsFrom = 'T';
		char bsTo = 'C';
		if (bsMapping && ((char *) extData)[i] == 1) {
			bsFrom = 'A';
			bsTo = 'G';
		}
		computeCigarMD(results[i], gpu_return_values[offset + alignments_per_thread * 3 + k], gpuCigar,
				refSeqList[i] + gpu_return_values[offset + k], qrySeqList[i], bsFrom, bsTo);
		results[i].PositionOffset = gpu_return_values[offset + k];
	}
	Log.Verbose("CPU Time: %.3fs", cpuTimer.ET());

	delete[] calignments;
	calignments = 0;
	delete[] gpu_return_values;
	gpu_return_values = 0;

#ifndef NDEBUG
	Log.Warning("Releasing results.");
#endif
	clReleaseMemObject(results_gpu);
#ifndef NDEBUG
	Log.Warning("Releasing alignments.");
#endif
	clReleaseMemObject(alignments_gpu);
#ifndef NDEBUG
	Log.Warning("Releasing matrix.");
#endif
	clReleaseMemObject(matrix_gpu);
	if (host->isGPU()) {
#ifndef NDEBUG
		Log.Warning("Releasing scaff.");
#endif
		clReleaseMemObject(c_scaff_gpu);
		if (bsMapping) {
			clReleaseMemObject(bsdirection_gpu);
		}
	}

//	//TODO: remove
//	for (int i = 0; i < batchSize; ++i) {
//		delete[] qrySeqList[i];
//	}
//	delete[] qrySeqList;
#ifndef NDEBUG
	Log.Message("SW finished computing alignments for %d sequences (elapsed: %.3fs)", batchSize_, timer.ET());
#endif

	delete[] tmpRefSeqList;
	delete[] tmpQrySeqList;

	return batchSize_;

}

int printCigarElement(char const op, short const length, char * cigar) {
	int offset = 0;
	offset = sprintf(cigar, "%d%c", length, op);
	return offset;
}

void debugCigar(int op, int length) {
	switch (op) {
	case CIGAR_M:
		Log.Message("CIGAR: %d M", length);
		break;
	case CIGAR_I:
		Log.Message("CIGAR: %d I", length);
		break;
	case CIGAR_D:
		Log.Message("CIGAR: %d D", length);
		break;
	case CIGAR_N:
		Log.Message("CIGAR: %d N", length);
		break;
	case CIGAR_S:
		Log.Message("CIGAR: %d S", length);
		break;
	case CIGAR_H:
		Log.Message("CIGAR: %d H", length);
		break;
	case CIGAR_P:
		Log.Message("CIGAR: %d P", length);
		break;
	case CIGAR_EQ:
		Log.Message("CIGAR: %d EQ", length);
		break;
	case CIGAR_X:
		Log.Message("CIGAR: %d X", length);
		break;
	default:
		Log.Error("Invalid cigar operator.");
		exit(1);
	}
}

void SWOclCigar::computeCigarMD(Align & result, int const gpuCigarOffset, short const * const gpuCigar, char const * const refSeq,
		char const * const qrySeq, char const bsFrom, char const bsTo) {

	static bool const bsMapping = Config.GetInt("bs_mapping") == 1;
	int cigar_offset = 0;
	int md_offset = 0;

	if ((gpuCigar[gpuCigarOffset] >> 4) > 0) {
		cigar_offset += printCigarElement('S', gpuCigar[gpuCigarOffset] >> 4, result.pCigar + cigar_offset);
		result.QStart = gpuCigar[gpuCigarOffset] >> 4;
	} else {
		result.QStart = 0;
	}

	int match = 0;
	int mismatch = 0;
	int total = 0;
	int cigar_m_length = 0;
	int md_eq_length = 0;
	int ref_index = 0;
	int read_index = result.QStart;
	for (int j = gpuCigarOffset + 1; j < (alignment_length - 1); ++j) {
		int op = gpuCigar[j] & 15;
		int length = gpuCigar[j] >> 4;

		//debugCigar(op, length);
		total += length;
		switch (op) {
		case CIGAR_X:
			cigar_m_length += length;
			if (!bsMapping)
				mismatch += length;

			//Produces: 	[0-9]+(([A-Z]+|\^[A-Z]+)[0-9]+)*
			//instead of: 	[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
			md_offset += sprintf(result.pMD + md_offset, "%d", md_eq_length);
			for (int k = 0; k < length; ++k) {
				if (bsMapping) {
					if (qrySeq[read_index] == bsFrom && refSeq[ref_index] == bsTo) {
						match += 1;
					} else {
						mismatch += 1;
					}
				}
				md_offset += sprintf(result.pMD + md_offset, "%c", refSeq[ref_index++]);
				read_index += 1;
			}
			md_eq_length = 0;

			break;
		case CIGAR_EQ:
			match += length;
			cigar_m_length += length;
			md_eq_length += length;
			ref_index += length;
			read_index += length;
			break;
		case CIGAR_D:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length, result.pCigar + cigar_offset);
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('D', length, result.pCigar + cigar_offset);

			md_offset += sprintf(result.pMD + md_offset, "%d", md_eq_length);
			md_eq_length = 0;
			result.pMD[md_offset++] = '^';
			for (int k = 0; k < length; ++k) {
				result.pMD[md_offset++] = refSeq[ref_index++];
			}

			break;
		case CIGAR_I:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length, result.pCigar + cigar_offset);
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('I', length, result.pCigar + cigar_offset);
			read_index += length;
			break;
		default:
			Log.Error("Invalid cigar string: %d", op);
			exit(1);
		}
	}
	md_offset += sprintf(result.pMD + md_offset, "%d", md_eq_length);
	if (cigar_m_length > 0) {
		cigar_offset += printCigarElement('M', cigar_m_length, result.pCigar + cigar_offset);
		cigar_m_length = 0;
	}

	if ((gpuCigar[alignment_length - 1] >> 4) > 0) {
		cigar_offset += printCigarElement('S', gpuCigar[alignment_length - 1] >> 4, result.pCigar + cigar_offset);
		result.QEnd = gpuCigar[alignment_length - 1] >> 4;
	} else {
		result.QEnd = 0;
	}

	result.pCigar[cigar_offset] = '\0';
	result.pMD[md_offset] = '\0';

	result.Identity = match * 1.0f / total;
	result.NM = mismatch;
}

long SWOclCigar::getMaxAllocSize(int const batch_size) {
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

	long alignments = (long) batch_size * (long) alignment_length * (long) 2 * (long) sizeof(cl_char);
	//	std::cout << "alignments: " << alignments << std::endl;
	return std::max(matrix, alignments);
}

int SWOclCigar::GetAlignBatchSize() const {
	return batch_size_align * step_count;
}

int SWOclCigar::computeAlignmentBatchSize() {

	cl_uint mpCount = host->getDeviceInfoInt(CL_DEVICE_MAX_COMPUTE_UNITS);
	//TODO: Fix
	//	unsigned long max_alloc = host->getDeviceInfoLong(CL_DEVICE_MAX_MEM_ALLOC_SIZE) * 2.0;
	int block_count = mpCount * block_multiplier * (host->getThreadPerMulti() / threads_per_block);
	block_count = (block_count / mpCount) * mpCount;

	unsigned long largest_alloc = getMaxAllocSize(block_count * threads_per_block);
	cl_mem testAlloc = 0;
	while (!host->testAllocate(largest_alloc)) {
		block_count -= mpCount;
		largest_alloc = getMaxAllocSize(block_count * threads_per_block);
#ifndef NDEBUG
		Log.Warning("Reducing batch size to %d", block_count * threads_per_block);
#endif
	}
	if (!host->isGPU()) {
		block_count *= 4;
	}
#ifndef NDEBUG
	Log.Message("Multi processor count: %d", mpCount);
	Log.Message("Max. threads per multi processor: %d", host->getThreadPerMulti());
	Log.Message("Threads per block used: %d", threads_per_block);
	Log.Message("Block number: %d", block_count);
	Log.Message("Batch size: %d", (block_count * threads_per_block));
	//TODO: Print debug info
#endif

	return block_count * threads_per_block;

}
