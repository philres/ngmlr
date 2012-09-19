/*
 * SWOcl.cpp
 *
 *  Created on: Apr 8, 2011
 *      Author: philipp_
 */

#include "SWOcl.h"
#include "OclHost.h"

#include <stdio.h>
#include <string.h>
#include <sstream>

#include "Timing.h"

//OpenCL source code
#include "oclDefines.h"
#include "oclSwScore.h"
#include "oclEndFreeScore.h"

using std::stringstream;

long seq_count = 0;
long overall = 0;

cl_program SWOcl::clProgram = 0;
int SWOcl::programUserCount = 0;

bool usedPinnedMemory = true;

int SWOcl::BatchScore(int const mode, int const batchSize_, char const * const * const refSeqList_, char const * const * const qrySeqList_,
		float * const results_, void * extData) {

	if (batchSize_ <= 0) {
		Log.Warning("Score for batchSize <= 0");
		return 0;
	}
#ifndef NDEBUG
	Log.Warning("BatchScore: %d", batchSize_);
	Log.Warning("Alignment mode: %d", mode);

	int count = 0;
	for (int i = 0; i < batchSize_; ++i) {
		//int count[256] = {};
		bool equal = true;
		char const * const ref = refSeqList_[i];
		for (int j = 0; j < config_ref_size - 1 && equal; ++j) {
			//count[ref[j]]+=1;
			equal = equal && ((ref[j] == ref[j + 1]) || ref[j + 1] == '\0');
		}
		if (equal) {
			count += 1;
			Log.Warning("Aligning read to %s", ref);
		}
	}
	Log.Error("Found %d reads homopolymer read.", count);
#endif

	bool batchSizeDif = !host->isGPU() && (batchSize_ % 4 != 0);
	int batchSize = (batchSizeDif) ? batchSize_ + 4 : batchSize_;
	float * const results = batchSizeDif ? new float[batchSize] : results_;
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

	cl_kernel scoreKernel = 0;
	switch (mode & 0xFF) {
	case 0: {
		scoreKernel = swScoringKernel;

		//Log.Verbose("Alignment mode: local");
	}
		break;
	case 1: {
		scoreKernel = swScoringKernelGlobal;

		//Log.Verbose("Alignment mode: end-free");
	}
		break;
	default:
		break;
		Log.Error("Unsupported alignment mode %i", mode & 0xFF);
		return 0;
	}
	Timer timer;
	timer.ST();

	cl_int ciErrNum = 0;
	cl_mem results_gpu = host->allocate(CL_MEM_WRITE_ONLY, batch_size_scoring * sizeof(cl_float));

	if (host->isGPU()) {
		ciErrNum |= clSetKernelArg(interleaveKernel, 0, sizeof(cl_mem), (void *) (&scaffold_gpu));
		ciErrNum |= clSetKernelArg(interleaveKernel, 1, sizeof(cl_mem), &c_scaff_gpu);

		ciErrNum |= clSetKernelArg(scoreKernel, 0, sizeof(cl_mem), (void *) (&c_scaff_gpu));
	} else {
		ciErrNum |= clSetKernelArg(scoreKernel, 0, sizeof(cl_mem), (void *) (&scaffold_gpu));
	}
	ciErrNum |= clSetKernelArg(scoreKernel, 1, sizeof(cl_mem), (void*) (&reads_gpu));
	ciErrNum |= clSetKernelArg(scoreKernel, 2, sizeof(cl_mem), &results_gpu);

	cl_mem bsdirection_gpu = 0;
	static bool const bsMapping = Config.GetInt("bs_mapping") == 1;
	if (bsMapping) {
		bsdirection_gpu = host->allocate(CL_MEM_READ_ONLY, batch_size_scoring * sizeof(cl_char));
		ciErrNum |= clSetKernelArg(scoreKernel, 3, sizeof(cl_mem), (void *) (&bsdirection_gpu));
	}

	runSwScoreKernel(scoreKernel, batchSize, qrySeqList, refSeqList, (char *) extData, bsdirection_gpu, results_gpu, results);

	clReleaseMemObject(results_gpu);
	if (bsMapping) {
		clReleaseMemObject(bsdirection_gpu);
	}
	seq_count += batchSize;

	overall += seq_count;
	seq_count = 0;

//#ifndef NDEBUG
	//Log.Verbose("SW finished computing score for %d sequences (elapsed: %.3fs)", batchSize, timer.ET());
//#endif

	if (batchSizeDif) {
		for (int i = 0; i < batchSize_; ++i) {
			results_[i] = results[i];
		}
		delete[] results;
	}
	delete[] tmpRefSeqList;
	delete[] tmpQrySeqList;

	return batchSize_;
}

SWOcl::SWOcl(char const * const oclSwScoreSourceCode, char const * const additional_defines, OclHost * phost) :
		host(phost) {

//	Log.Message("HHHHHHHHHHHHHAAAALLLOOO");
//	if (!host->checkLocalMemory(46000)) {
//		Log.Warning("Reducing threads per block to 128 due to lack of local memory.");
//		threads_per_block = 256;
//	} else {
//		threads_per_block = 256;
//	}

	alignments_per_thread = 1;
	if (!host->isGPU()) {
		alignments_per_thread = 4;
	}

	//block_count = 0;
	batch_size_scoring = 0;
	//batch_size_align = 0;
	alignment_length = 0;
	matrix_length = 0;
	config_ref_size = -1;

	read_data_size = 0;
	cpu_read_data = 0;
	ref_data_size = 0;
	cpu_ref_data = 0;

	Timer timer2;
	timer2.ST();
#ifndef NDEBUG
	Log.Message("SW finished init opencl. (elapsed: %.3fs)", timer2.ET());
#endif

	Timer timer;
	timer.ST();

	checkMemory();

	//Log.Message("Source: %s", oclSwScoreSourceCode);

	stringstream buildCmd;
	static bool const bsMapping = Config.GetInt("bs_mapping") == 1;
	buildCmd << "-D MATRIX_LENGTH=" << (matrix_length * threads_per_block) << " -D interleave_number=" << "256" << " -D threads_per_block="
			<< threads_per_block << " -D match=" << Config.GetFloat("score_match") << " -D mismatch=" << Config.GetFloat("score_mismatch")
			<< " -D gap_read=" << Config.GetFloat("score_gap_read") << " -D gap_ref=" << Config.GetFloat("score_gap_ref") << " -D matchBS="
			<< Config.GetFloat("score_match_tt") << " -D mismatchBS=" << Config.GetFloat("score_mismatch_tc") << " -D read_length="
			<< Config.GetInt("qry_max_len") << " -D ref_length=" << config_ref_size << " -D corridor_length="
			<< (Config.GetInt("corridor") + 1) << " -D alignment_length=" << alignment_length << " " << additional_defines;
	if (host->isGPU()) {
		buildCmd << " -D __GPU__";
	} else {
		buildCmd << " -D __CPU__";
	}
	if (bsMapping) {
		buildCmd << " -D __BS__";
	}

#pragma omp critical
	{
		if (programUserCount == 0) {
#ifndef NDEBUG
			Log.Message("Building program.");
#endif
			stringstream source;
			source << oclDefines << std::endl << oclSwScore << oclEndFreeScore << oclSwScoreSourceCode;
			//			Log.Message("SOURCE: %s", oclEndFreeScore);
			clProgram = host->setUpProgram(source.str().c_str(), buildCmd.str().c_str());
		}
		programUserCount += 1;
	}
	swScoringKernel = host->setupKernel(clProgram, "oclSW");
	swScoringKernelGlobal = host->setupKernel(clProgram, "oclSW_Global");
	interleaveKernel = host->setupKernel(clProgram, "interleaveSeq");
	//Log.Message("Releasing program");

	if (host->isGPU()) {

		reads_gpu = host->allocate(CL_MEM_READ_ONLY, read_data_size);
		scaffold_gpu = host->allocate(CL_MEM_READ_WRITE, ref_data_size);
		c_scaff_gpu = host->allocate(CL_MEM_READ_WRITE, ref_data_size * sizeof(cl_char));

		if (usedPinnedMemory) {
			//Setup pinned memory
			cmPinnedBufRead = host->allocate(CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, read_data_size);
			cmPinnedBufRef = host->allocate(CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, ref_data_size);
			cpu_read_data = (char*) host->mapBuffer(cmPinnedBufRead, 0, read_data_size);
			cpu_ref_data = (char*) host->mapBuffer(cmPinnedBufRef, 0, ref_data_size);
		} else {
			cpu_read_data = new char[read_data_size];
			cpu_ref_data = new char[ref_data_size];
		}
	} else {
		cpu_read_data = new char[read_data_size];
		cpu_ref_data = new char[ref_data_size];
		reads_gpu = host->allocate(CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, read_data_size, cpu_read_data);
		scaffold_gpu = host->allocate(CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, ref_data_size, cpu_ref_data);
	}
#ifndef NDEBUG
	Log.Message("SW finished allocating memory. (elapsed: %.3fs)", timer.ET());
#endif
}

SWOcl::~SWOcl() {

	//#ifndef NDEBUG
	//	Log.Warning("Overall score count %d", overall);
	//	Timer timer;
	//	timer.ST();
	//#endif
	//
	//
	clReleaseKernel(swScoringKernel);
	clReleaseKernel(swScoringKernelGlobal);
	clReleaseKernel(interleaveKernel);
	if (--programUserCount == 0) {
#ifndef NDEBUG
		Log.Message("Releasing ocl program.");
#endif
		clReleaseProgram(clProgram);
		clProgram = 0;
	}
	clReleaseMemObject(reads_gpu);
	clReleaseMemObject(scaffold_gpu);
	if (host->isGPU()) {
		clReleaseMemObject(c_scaff_gpu);
		if (usedPinnedMemory) {
			Log.Message("Releasing pinned memory.");
			clReleaseMemObject(cmPinnedBufRead);
			clReleaseMemObject(cmPinnedBufRef);
		}
	} else {
		if (cpu_read_data != 0) {
			delete[] cpu_read_data;
			cpu_read_data = 0;
		}
		if (cpu_ref_data != 0) {
			delete[] cpu_ref_data;
			cpu_ref_data = 0;
		}
	}
	//
	//#ifndef NDEBUG
	//	Log.Message("SW finished releasing memory. (elapsed: %.3fs)", timer.ET());
	//#endif
}

int SWOcl::getMaxAllocSize() {
	return std::max(std::max(ref_data_size * sizeof(cl_char), read_data_size * sizeof(cl_char)), batch_size_scoring * sizeof(cl_float));
}

void SWOcl::checkMemory() {

	//Get Paramter
	block_multiplier = std::max(Config.GetInt("block_multiplier"), 1);
	/*	if(!host->isGPU()) {
	 block_multiplier = 1;
	 }*/

#ifdef _TEST
	step_count = 1;
#else
	step_count = std::max(Config.GetInt("step_count"), 1);
	if (!host->isGPU()) {
		step_count = 1;
	}
#endif

	alignment_length = 2 * Config.GetInt("qry_max_len") + Config.GetInt("corridor") + 1;
	config_ref_size = Config.GetInt("qry_max_len") + Config.GetInt("corridor");
	matrix_length = Config.GetInt("corridor") + 1;
	batch_size_scoring = computeScoringBatchSize();
	read_data_size = batch_size_scoring * Config.GetInt("qry_max_len");
	ref_data_size = batch_size_scoring * config_ref_size;

	//Check local Memory
	unsigned int localMemByte = (host->isGPU()) ? (matrix_length * threads_per_block) * 2 : 0;
	//cl_ulong localMemAvailable = host->getDeviceInfoLong(CL_DEVICE_LOCAL_MEM_SIZE);
#ifndef NDEBUG
	Log.Message("You are using %d byte local memory. (%d byte available).", localMemByte, host->getDeviceInfoLong(CL_DEVICE_LOCAL_MEM_SIZE));
#endif
	if (!host->checkLocalMemory(localMemByte)) {
		Log.Error("Not enough local memory available. Please reduce corridor size.");
		Log.Message("You are using %d byte local memory. (%d byte available).", localMemByte, host->getDeviceInfoLong(CL_DEVICE_LOCAL_MEM_SIZE));
		exit(1);
	}

	//Check global Memory
	//unsigned int gMemScoreGPU = read_data_size + ref_data_size + batch_size_scoring * sizeof(cl_float);
	//unsigned int gMemScoreCPU = read_data_size + ref_data_size;

	//unsigned int gMemAlignGPU = read_data_size + ref_data_size + batch_size_align * sizeof(cl_float) + 3 * batch_size_align * sizeof(cl_short) + batch_size_align * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_short) + batch_size_align * alignment_length * 2 * sizeof(cl_short);
	//unsigned int gMemAlignCPU = read_data_size + ref_data_size + batch_size_align * alignment_length * 2 * sizeof(short) + 3 * batch_size_align * sizeof(short);

	//cl_ulong gMemAvailable = host->getDeviceInfoLong(CL_DEVICE_GLOBAL_MEM_SIZE) * 0.9;

	//#ifndef NDEBUG
	//Log.Message("Scoring:");
	//Log.Message("For a batch size of %d, %d byte global memory will be required. The selected device has %d bytes available.", batch_size_scoring, gMemScoreGPU, gMemAvailable);
	//Log.Message("Additional RAM needed: %d byte", gMemScoreCPU);
	//Log.Message("Alignment:");
	//Log.Message("For a batch size of %d, %d byte global memory will be required. The selected device has %d bytes available.", batch_size_align, gMemAlignGPU, gMemAvailable);
	//Log.Message("Additional RAM needed: %d byte", gMemAlignCPU);
	//#endif
	//	if (gMemAlignGPU > gMemAvailable || gMemScoreGPU > gMemAvailable) {
	//		Log.Error("Not enough global memory available. Please reduce batch size.");
	//		exit(1);
	//	}
}

int SWOcl::GetScoreBatchSize() const {
	return batch_size_scoring * step_count;
}

//void flatten(char const * const * const src, char * dest, int const size, int const alignment_number, int const batch_size) {
//
//	for (int i = 0; i < alignment_number; ++i) {
//		memcpy64(&dest[i * size], src[i], sizeof(char) * size);
//	}
//	if (alignment_number < batch_size) { // if there are not enough sequence pairs
//		//set the rest
//		//Log.Warning("Batchsize not large enough!");
//		memset(&dest[alignment_number * size], '\0', (batch_size - alignment_number) * size * sizeof(char));
//	}
//}

#ifdef _TEST
#include "Timer.h"
#define benchCPP(x,z)  test.start(); for(int y=0;y<100;++y) x; test.stop(); Log.Message("%s\t%f", z, test.getElapsedTimeInMilliSec());
#define benchOCL(x,z)  test.start(); for(int y=0;y<100;++y) x;host->waitForDevice(); test.stop(); Log.Message("%s\t%f", z, test.getElapsedTimeInMilliSec());
#else
#define benchCPP(x,z) x;
#define benchOCL(x,z) x;
#endif

int SWOcl::runSwScoreKernel(cl_kernel scoreKernel, const int batchSize, const char * const * const qrySeqList,
		const char * const * const refSeqList, char * bsDirection, cl_mem & bsdirection_gpu, cl_mem & results_gpu, float * const results) {
#ifdef _TEST
	Timer2 test;
#endif

	const size_t cnDim = batch_size_scoring / alignments_per_thread;
	const size_t cBlockSize = threads_per_block;
	int runbatchSize = std::min(batch_size_scoring, batchSize);

	benchCPP( copySeqDataToDevice(cpu_read_data, cpu_ref_data, qrySeqList, refSeqList, runbatchSize, batch_size_scoring), "Flatten: ");

	for (int i = 0; i < batchSize; i += batch_size_scoring) {
		if (host->isGPU()) {
			benchOCL( host->writeToDevice(scaffold_gpu, CL_FALSE, 0, ref_data_size * sizeof(cl_char), cpu_ref_data), "Scaff: ");
			benchOCL( host->writeToDevice(reads_gpu, CL_FALSE, 0, read_data_size * sizeof(cl_char), cpu_read_data), "Read: ");
			if (bsdirection_gpu != 0) {
				host->writeToDevice(bsdirection_gpu, CL_FALSE, 0, std::min(runbatchSize, batchSize - i) * sizeof(cl_char), bsDirection + i);
			}
			host->waitForDevice();
			benchOCL(host->executeKernel(interleaveKernel, cnDim, cBlockSize), "Interl: ");
		}
		benchOCL(host->executeKernel(scoreKernel, cnDim, cBlockSize);, "Kernel: ");

		int nextBatch = (i + batch_size_scoring);
		int nextRunbatchSize = std::min(batch_size_scoring, batchSize - nextBatch);
		if (nextRunbatchSize > 0) {
			copySeqDataToDevice(cpu_read_data, cpu_ref_data, qrySeqList + nextBatch, refSeqList + nextBatch, nextRunbatchSize,
					batch_size_scoring);
#ifdef _TEST
			Log.Error("I should not be here. Set step_count to 1 when compiling with _TEST switch.");
			throw "";
#endif
		}

		benchOCL( host->readFromDevice(results_gpu, CL_FALSE, 0, runbatchSize, results + i, sizeof(cl_float)), "Result: ");
		runbatchSize = nextRunbatchSize;
	}
	host->waitForDevice();

	return 0;
}

//void SWOcl::copy(char const * const * const refs, char * cpu_ref_data, int const alignment_number, int const batch_size, int const refsize) {
//
//	int const size = 16;
//	int offset = -1;
//
//	for (int j = 0; j < alignment_number; j += size) {
//		for (int readPos = 0; readPos < refsize; ++readPos) {
//			for (int i = 0; i < size; ++i) {
//				cpu_ref_data[++offset] = ((i + j) < alignment_number) ? refs[j
//						+ i][readPos] : '\0';
//			}
//		}
//	}
//}

//void SWOcl::copy(char const * const * const refs, char * cpu_ref_data, int const alignment_number, int const batch_size, int const refsize) {
//
//	int const size = 16;
//	int offset = -1;
//
//	for (int j = 0; j < alignment_number; j += size) {
//		for (int readPos = 0; readPos < refsize; ++readPos) {
//			for (int i = 0; i < size; ++i) {
//				cpu_ref_data[++offset] = refs[j + i][readPos];
//			}
//		}
//	}
//	//memset(cpu_ref_data + offset, '\0', (refsize * batch_size) - offset);
//}

void SWOcl::copy(char const * const * const refs, char * cpu_ref_data, int const alignment_number, int const batch_size,
		int const refsize) {

	int const size = 16;
	int offset = 0;

	int j = 0;
	for (j = 0; j < (alignment_number / size) * size; j += size) {
		for (int readPos = 0; readPos < refsize; ++readPos) {
			cpu_ref_data[offset] = refs[j][readPos];
			cpu_ref_data[offset + 1] = refs[j + 1][readPos];
			cpu_ref_data[offset + 2] = refs[j + 2][readPos];
			cpu_ref_data[offset + 3] = refs[j + 3][readPos];
			cpu_ref_data[offset + 4] = refs[j + 4][readPos];
			cpu_ref_data[offset + 5] = refs[j + 5][readPos];
			cpu_ref_data[offset + 6] = refs[j + 6][readPos];
			cpu_ref_data[offset + 7] = refs[j + 7][readPos];
			cpu_ref_data[offset + 8] = refs[j + 8][readPos];
			cpu_ref_data[offset + 9] = refs[j + 9][readPos];
			cpu_ref_data[offset + 10] = refs[j + 10][readPos];
			cpu_ref_data[offset + 11] = refs[j + 11][readPos];
			cpu_ref_data[offset + 12] = refs[j + 12][readPos];
			cpu_ref_data[offset + 13] = refs[j + 13][readPos];
			cpu_ref_data[offset + 14] = refs[j + 14][readPos];
			cpu_ref_data[offset + 15] = refs[j + 15][readPos];
			offset += 16;
		}
	}
	if (j < alignment_number) {
		for (int readPos = 0; readPos < refsize; ++readPos) {
			cpu_ref_data[offset] = (j < alignment_number) ? refs[j][readPos] : '\0';
			cpu_ref_data[offset + 1] = ((j + 1) < alignment_number) ? refs[j + 1][readPos] : '\0';
			cpu_ref_data[offset + 2] = ((j + 2) < alignment_number) ? refs[j + 2][readPos] : '\0';
			cpu_ref_data[offset + 3] = ((j + 3) < alignment_number) ? refs[j + 3][readPos] : '\0';
			cpu_ref_data[offset + 4] = ((j + 4) < alignment_number) ? refs[j + 4][readPos] : '\0';
			cpu_ref_data[offset + 5] = ((j + 5) < alignment_number) ? refs[j + 5][readPos] : '\0';
			cpu_ref_data[offset + 6] = ((j + 6) < alignment_number) ? refs[j + 6][readPos] : '\0';
			cpu_ref_data[offset + 7] = ((j + 7) < alignment_number) ? refs[j + 7][readPos] : '\0';
			cpu_ref_data[offset + 8] = ((j + 8) < alignment_number) ? refs[j + 8][readPos] : '\0';
			cpu_ref_data[offset + 9] = ((j + 9) < alignment_number) ? refs[j + 9][readPos] : '\0';
			cpu_ref_data[offset + 10] = ((j + 10) < alignment_number) ? refs[j + 10][readPos] : '\0';
			cpu_ref_data[offset + 11] = ((j + 11) < alignment_number) ? refs[j + 11][readPos] : '\0';
			cpu_ref_data[offset + 12] = ((j + 12) < alignment_number) ? refs[j + 12][readPos] : '\0';
			cpu_ref_data[offset + 13] = ((j + 13) < alignment_number) ? refs[j + 13][readPos] : '\0';
			cpu_ref_data[offset + 14] = ((j + 14) < alignment_number) ? refs[j + 14][readPos] : '\0';
			cpu_ref_data[offset + 15] = ((j + 15) < alignment_number) ? refs[j + 15][readPos] : '\0';
			offset += 16;
		}
	}
}

void SWOcl::copySeqDataToDevice(char * cpu_read_data, char * cpu_ref_data, char const * const * const reads,
		char const * const * const refs, int const alignment_number, int const batch_size) {
	int refsize = Config.GetInt("qry_max_len") + Config.GetInt("corridor");

	int readsize = Config.GetInt("qry_max_len");

	//memset(cpu_ref_data, '\0', ref_data_size);
	//memset(cpu_read_data, '\0', read_data_size);

	//copy(refs, cpu_ref_data, alignment_number, batch_size, refsize);
	//copy(reads, cpu_read_data, alignment_number, batch_size, readsize);

	for (int i = 0; i < alignment_number; ++i) {
		memcpy(&cpu_ref_data[i * config_ref_size], refs[i], sizeof(char) * refsize);
		memcpy(&cpu_read_data[i * readsize], reads[i], sizeof(char) * readsize);
	}
	if (alignment_number < batch_size) { // if there are not enough sequence pairs
		//set the rest
		//Log.Warning("Batchsize not large enough!");
		memset(&cpu_read_data[alignment_number * readsize], '\0', (batch_size - alignment_number) * readsize * sizeof(char));
		memset(&cpu_ref_data[alignment_number * config_ref_size], '\0', (batch_size - alignment_number) * config_ref_size * sizeof(char));
	}

}

unsigned int SWOcl::computeScoringBatchSize() {

	cl_uint mpCount = host->getDeviceInfoInt(CL_DEVICE_MAX_COMPUTE_UNITS);
	//long max_alloc = host->getDeviceInfoLong(CL_DEVICE_MAX_MEM_ALLOC_SIZE);
	int block_count = mpCount * block_multiplier * (host->getThreadPerMulti() / threads_per_block);
	block_count = (block_count / mpCount) * mpCount;

	//long largest_alloc = (block_count * threads_per_block) * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char);
	//while (largest_alloc > max_alloc) {
	//		Log.Warning("Reducing batch size.");
	//block_count -= mpCount;
	//largest_alloc = (block_count * threads_per_block) * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char);
	//}
#ifndef NDEBUG
	Log.Message("Multi processor count: %d", mpCount);
	Log.Message("Max. threads per multi processor: %d", host->getThreadPerMulti());
	Log.Message("Threads per block used: %d", threads_per_block);
	Log.Message("Block number: %d", block_count);
	Log.Message("Batch size: %d", (block_count * threads_per_block * alignments_per_thread));
	//TODO: Print debug info
#endif

	return block_count * threads_per_block * alignments_per_thread;
}

