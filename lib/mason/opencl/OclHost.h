/*
 * OclHost.h
 *
 *  Created on: May 25, 2011
 *      Author: philipp_
 */

#ifndef OCLHOST_H_
#define OCLHOST_H_

//#ifdef __APPLE__
//#include <OpenCL/opencl.h>
//#else
#include <CL/opencl.h>
#//endif

#include <string>

#include "ILog.h"

class OclHost {

	private:
		cl_device_type const devType;
		cl_context partitionDevice(cl_platform_id platform, cl_uint ciDeviceCount, cl_device_id *cdDevices, cl_int cores);
//		cl_device_id getDevice(cl_context context, unsigned int gpu_id);
    cl_platform_id getPlatform();
	public:
		OclHost();
		OclHost(int const device_type, int gpu_id, int const cpu_cores);
		virtual ~OclHost();
		void initOpenCL(unsigned int gpu_id);
//		cl_device_id getDevice(unsigned int gpu_id);
		char *print_cl_errstring(cl_int err);
		const char *print_cl_buildstatus(cl_build_status s);
		void checkClError(const char *msg, cl_int ciErrNum);
		cl_program setUpProgram(const char * const oclSwScore, std::string buildOptions);
		cl_kernel setupKernel(cl_program program, const char * const kernelName);

		virtual int getThreadPerMulti();

		cl_mem allocate(cl_mem_flags, size_t size, void * ptr = NULL);

		bool checkGlobalMemory(size_t const size);

		bool checkLocalMemory(size_t const size);

		bool testAllocate(unsigned long size);

		void * mapBuffer(cl_mem buffer, size_t offset, size_t size);

		void writeToDevice(cl_mem buffer, cl_bool blocking_write, size_t offset, size_t size, const void * ptr);

		void waitForDevice();

		void executeKernel(cl_kernel kernel, const size_t global_work_size, const size_t local_work_size);

		void readFromDevice(cl_mem buffer, cl_bool blocking_read, size_t offset, size_t size, void * ptr, size_t size_of);

		bool isGPU();

		cl_uint getDeviceInfoInt(cl_device_info info);
		cl_ulong getDeviceInfoLong(cl_device_info info);

		static int contextUserCount;
		static cl_context oclGpuContext;
		static cl_device_id * devices;
		static cl_uint ciDeviceCount;

		//cl_device_id * devices;
//		cl_context oclGpuContext;
		cl_command_queue oclCommandQueue;
		cl_device_id oclDevice;

		cl_ulong maxGlobalMem;
		cl_ulong maxLocalMem;

};

#endif /* OCLHOST_H_ */
