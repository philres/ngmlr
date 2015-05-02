/*
 * OclHost.cpp
 *
 *  Created on: May 25, 2011
 *      Author: philipp_
 */

#include "OclHost.h"

#include <string.h>
#include <pthread.h>
#include "IConfig.h"

pthread_mutex_t mutext_next_sub_block;

#undef module_name
#define module_name "OPENCL"

char clPlatformName[1024];

//static clCreateSubDevicesEXT_fn pfn_clCreateSubDevicesEXT = NULL;

cl_context OclHost::oclGpuContext = 0;
int OclHost::contextUserCount = 0;
cl_device_id * OclHost::devices = 0;
cl_uint OclHost::ciDeviceCount = 0;

enum PLATFORM {
	NVIDIA = 0, AMD = 1
};

PLATFORM platform;

OclHost::OclHost() :
		devType(CL_DEVICE_TYPE_GPU) {
	OclHost(0, 0, 1);

}

OclHost::OclHost(int const device_type, int gpu_id, int const cpu_cores) :
		devType(device_type), maxGlobalMem(0), maxLocalMem(0) {
//		if (!isGPU()) {
//				gpu_id = 0;
//		}

	cl_int ciErrNum = CL_SUCCESS;
	Log.Verbose("Using device number %d", gpu_id);
//#pragma omp critical
//	{
	if (contextUserCount == 0) {
		Log.Verbose("Creating ocl context.");
//		cl_uint ciDeviceCount = 0;
		cl_platform_id cpPlatform = NULL;

		cpPlatform = getPlatform();
		//Get the devices

		//Get number of devices
		ciErrNum = clGetDeviceIDs(cpPlatform, devType, 0, NULL, &ciDeviceCount);
		checkClError("Couldn't get number of OpenCl devices. Error: ",
				ciErrNum);

		if (isGPU()) {
			//Getting device ids
			devices = (cl_device_id *) malloc(
					ciDeviceCount * sizeof(cl_device_id));
			ciErrNum = clGetDeviceIDs(cpPlatform, devType, ciDeviceCount,
					devices, NULL);
			checkClError("Couldn't get OpenCl device ids. Error: ", ciErrNum);

			//Create context
			oclGpuContext = clCreateContext(0, ciDeviceCount, devices, NULL,
					NULL, &ciErrNum);
			checkClError("Couldn't create context. Error: ", ciErrNum);
			Log.Message("Context for GPU devices created.");

			Log.Message("%d GPU device(s) found: ", ciDeviceCount);
			for (int i = 0; i < ciDeviceCount; ++i) {
				char device_string[1024];
				char driver_string[1024];
				clGetDeviceInfo(devices[i], CL_DEVICE_NAME,
						sizeof(device_string), &device_string, NULL);
				clGetDeviceInfo(devices[i], CL_DRIVER_VERSION,
						sizeof(driver_string), &driver_string, NULL);
				Log.Message("Device %d: %s (Driver: %s)", i, device_string, driver_string);
			}

		} else {
			if (ciDeviceCount > 1) {
				Log.Error("More than one CPU device found.");
				exit(-1);
			}

			cl_device_id device_id;
			ciErrNum = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_CPU, 1,
					&device_id, NULL);
			checkClError("Couldn't get CPU device id. Error: ", ciErrNum);

			Log.Message("%d CPU device found.", ciDeviceCount);
			char device_string[1024];
			char driver_string[1024];
			clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(device_string),
					&device_string, NULL);
			clGetDeviceInfo(device_id, CL_DRIVER_VERSION, sizeof(driver_string),
					&driver_string, NULL);
			Log.Message("Device %d: %s (Driver: %s)", 0, device_string, driver_string);

			cl_device_partition_property props[3];

			props[0] = CL_DEVICE_PARTITION_EQUALLY; // Equally
			props[1] = 1; // 4 compute units per sub-device
			props[2] = 0;

			devices = (cl_device_id *) malloc(256 * sizeof(cl_device_id));
			ciErrNum = clCreateSubDevices(device_id, props, 256, devices,
					&ciDeviceCount);
			if (ciErrNum == -18) {
				ciDeviceCount = 1;
				devices[0] = device_id;
			} else {
				checkClError("Couldn't create sub-devices. Error: ", ciErrNum);
			}

			Log.Message("%d CPU cores available.", ciDeviceCount);

			//Create context
			oclGpuContext = clCreateContext(0, ciDeviceCount, devices, NULL,
					NULL, &ciErrNum);
			checkClError("Couldn't create context. Error: ", ciErrNum);

		}
	}
	contextUserCount += 1;
	//}

	if (!isGPU()) {
		gpu_id = gpu_id % ciDeviceCount;
	}
	oclDevice = devices[gpu_id];
	//Create context
	//oclGpuContext = clCreateContext(0, 1, &oclDevice, NULL, NULL, &ciErrNum);
	//checkClError("Couldn't create context. Error: ", ciErrNum);

	// create command queue
	oclCommandQueue = clCreateCommandQueue(oclGpuContext, oclDevice, 0,
			&ciErrNum);

	checkClError("Couldn't create command queue for device: ", ciErrNum);

}

OclHost::~OclHost() {
	clReleaseCommandQueue(oclCommandQueue);
	//clReleaseDeviceEXT(oclDevice);

	if (--contextUserCount == 0) {
		Log.Verbose("Releasing ocl context.");
		clReleaseContext(oclGpuContext);
		oclGpuContext = 0;
	}

}

void OclHost::checkClError(char const * msg, cl_int ciErrNum) {
	if (ciErrNum != CL_SUCCESS) {
		Log.Error("%s\nError: %s (%d)", msg, print_cl_errstring(ciErrNum), ciErrNum);
		throw;
	}
}

cl_platform_id getPlatformID(char const * const platformName) {
	cl_int ciErrNum = 0;

	cl_uint platformNumber = 0;
	ciErrNum |= clGetPlatformIDs(0, 0, &platformNumber);
	if (ciErrNum == CL_SUCCESS) {

		//cdDevices = (cl_device_id *) malloc(ciDeviceCount
		//				* sizeof(cl_device_id));
		// Get OpenCL platform count
		cl_platform_id clPlatformID[platformNumber];
		ciErrNum = clGetPlatformIDs(platformNumber, clPlatformID, 0);
		Log.Message("Available platforms: %d", platformNumber);
		for (size_t i = 0; i < platformNumber; ++i) {
			ciErrNum = clGetPlatformInfo(clPlatformID[i], CL_PLATFORM_NAME,
					1024, &clPlatformName, NULL);
			Log.Message("%s", clPlatformName);
			if (ciErrNum == CL_SUCCESS) {
				if (strcasestr(clPlatformName, platformName) != 0) {
					Log.Message("Selecting OpenCl platform: %s", clPlatformName);
					ciErrNum = clGetPlatformInfo(clPlatformID[i],
					CL_PLATFORM_VERSION, 1024, &clPlatformName, NULL);
					Log.Message("Platform: %s", clPlatformName);
					//ciErrNum = clGetDeviceInfo(oclDevice, CL_DRIVER_VERSION, 1024, &clPlatformName, NULL);
					//ciErrNum = clGetPlatformInfo(clPlatformID[i], CL_DRIVER_VERSION, 1024, &clPlatformName, NULL);
					//Log.Message("Driver: %s", clPlatformName);
					return clPlatformID[i];
				}
			} else {
				Log.Error("Couldn't get OpenCl platform name. Error: ", ciErrNum);
			}
		}
	} else {
		Log.Error("Couldn't get OpenCl platform ids. Error: ", ciErrNum);
	}
	return 0;
}

cl_context OclHost::partitionDevice(cl_platform_id platform,
		cl_uint ciDeviceCount, cl_device_id *cdDevices, cl_int cores) {
	cl_uint numSubDevices = 0;
	cl_int ciErrNum = 0;
//	cl_context oclCPUContext = clCreateContext(0, ciDeviceCount, cdDevices,
//			NULL, NULL, &ciErrNum);
	//cl_device_id device_id = 0;

	Log.Message("%d", ciDeviceCount);
//	clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device_id, &ciDeviceCount);
//	Log.Message("%d", ciDeviceCount);

	//pfn_clCreateSubDevicesEXT = (clCreateSubDevicesEXT_fn) (clGetExtensionFunctionAddress("clCreateSubDevicesEXT"));
	cl_device_partition_property partitionPrty[3];

	partitionPrty[0] = CL_DEVICE_PARTITION_EQUALLY;
	partitionPrty[1] = 1;
	partitionPrty[2] = 0;

//	pfn_clCreateSubDevicesEXT(cdDevices[0], partitionPrty, 0, NULL, &numSubDevices);
	clCreateSubDevices(cdDevices[0], partitionPrty, 0, NULL, &numSubDevices);
	Log.Message("%d", numSubDevices);
	cl_device_id *subDevices = (cl_device_id*) (malloc(
			numSubDevices * sizeof(cl_device_id)));
	clCreateSubDevices(cdDevices[0], partitionPrty, numSubDevices, subDevices,
			NULL);
	// Create context for sub-devices
	cl_context context = clCreateContext(0, 1, subDevices, NULL, NULL,
			&ciErrNum);
	checkClError("BLABLABLAB", ciErrNum);
	Log.Verbose("Dividing CPU into %d devices.", numSubDevices);
	free(subDevices);
//	clReleaseDevice(device);
//	clReleaseContext(oclCPUContext);
	return context;
}

cl_platform_id OclHost::getPlatform() {
	cl_platform_id cpPlatform;
#ifdef __APPLE__
	cpPlatform = getPlatformID("APPLE");
	if (cpPlatform == 0) {
		Log.Error("No OpenCl platform found.");
		exit(1);
	}
#else
	//Get the first platform
	if (isGPU()) {
		cpPlatform = getPlatformID("NVIDIA");
		platform = NVIDIA;
		if (cpPlatform == 0) {
			Log.Warning("NVIDIA Platform not found. Falling back to AMD.");
			cpPlatform = getPlatformID("AMD");
			platform = AMD;
			if (cpPlatform == 0) {
				Log.Error("No OpenCl platform found.");
				exit(1);
			}
		}
	} else {
		cpPlatform = getPlatformID("AMD");
		if (cpPlatform == 0) {
			Log.Warning("AMD Platform not found. Falling back to Intel.");
			cpPlatform = getPlatformID("Intel");
			if (cpPlatform == 0) {
				Log.Error("No OpenCl platform found.");
				exit(1);
			}
		}
	}
#endif
	return cpPlatform;
}

void OclHost::initOpenCL(unsigned int cores) {
	//pthread_mutex_lock(&mutext_next_sub_block);

	//pthread_mutex_unlock(&mutext_next_sub_block);
}

bool OclHost::checkGlobalMemory(size_t const size) {
	return true;
	return (size
			< ((maxGlobalMem != 0) ? maxGlobalMem : (maxGlobalMem =
												getDeviceInfoLong(
												CL_DEVICE_MAX_MEM_ALLOC_SIZE))));
}

bool OclHost::checkLocalMemory(size_t const size) {
	return (size
			< ((maxLocalMem != 0) ?
					maxLocalMem :
					(maxLocalMem = getDeviceInfoLong(CL_DEVICE_LOCAL_MEM_SIZE))));
}

#include <iostream>

bool OclHost::testAllocate(unsigned long size) {
	//cl_ulong const maxBlockSize = (getDeviceInfoLong(CL_DEVICE_GLOBAL_MEM_SIZE) / 1024 * 1000) / (unsigned long)(Config.GetInt("cpu_threads") * 8);
	//std::cout << "MaxBlockSize: " << maxBlockSize << std::endl;
	//if(size < maxBlockSize) {
	cl_ulong gMem = (getDeviceInfoLong(CL_DEVICE_GLOBAL_MEM_SIZE)) * 0.45;
	if (platform == AMD) {
		return size < (getDeviceInfoLong(CL_DEVICE_MAX_MEM_ALLOC_SIZE))
				&& (size * (unsigned long) Config.GetInt("cpu_threads")) < gMem;
	}

	if ((size * (unsigned long) Config.GetInt("cpu_threads")) < gMem) {
		cl_int errCode = 0;

		//return true;
		cl_mem mem = clCreateBuffer(oclGpuContext, 0, size, 0, &errCode);
		if (mem != 0 && errCode == CL_SUCCESS) {
			clReleaseMemObject(mem);
			//		delete[] test;
			//		test = 0;
			return true;
		} else {
			return false;
		}
	}
	//}
	return false;
}

cl_mem OclHost::allocate(cl_mem_flags flags, size_t size, void * ptr) {
	cl_int errCode = 0;
	cl_mem mem = 0;
	if (checkGlobalMemory(size)) {
		Log.Verbose("Allocationg %d bytes on opencl device.", size);
		mem = clCreateBuffer(oclGpuContext, flags, size, ptr, &errCode);
		checkClError("Unable to create buffer.", errCode);
	} else {
		Log.Error("Unable to allocate %d byte of global memory. Max allocation size is %d on your device. Please reduce batch size.", size, maxGlobalMem);
		return 0;
	}
	return mem;
}

void * OclHost::mapBuffer(cl_mem buffer, size_t offset, size_t size) {
	cl_int errCode = 0;
	//Log.Verbose("Pinning %d memory.", size);
	void * ptr = clEnqueueMapBuffer(oclCommandQueue, buffer, CL_TRUE,
	CL_MAP_WRITE, offset, size, 0, NULL, NULL, &errCode);
	//clFlush(oclCommandQueue);
	checkClError("Unable to map buffer.", errCode);
	return ptr;
}

cl_kernel OclHost::setupKernel(cl_program program,
		char const * const kernelName) {
	// Create Kernel
	cl_int ciErrNum = 0;
	cl_kernel kernel = clCreateKernel(program, kernelName, &ciErrNum);

	//size_t test = 0;
	//clGetKernelWorkGroupInfo(kernel, oclDevice, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &test,0);
	//Log.Message("Kernel workgroup size multiple: %d", test);

	checkClError("Unable to create OpenCl kernel.", ciErrNum);
	return kernel;
}

cl_program OclHost::setUpProgram(char const * const oclSwScore,
		std::string buildOptions) {

	//char const * additional_options_nv = " -cl-nv-verbose -cl-fast-relaxed-math";
	//char * additional_options = 0;

	size_t program_length = strlen(oclSwScore);
	if (strcasestr(clPlatformName, "NVIDIA") != 0) {
		buildOptions += " -D __NVIDIA__ -cl-nv-verbose -cl-fast-relaxed-math";
	}
	cl_int ciErrNum = 0;
	//Log.Message("Source: %s\n===========================", oclSwScore);
	// create the program
	cl_program cpProgram = clCreateProgramWithSource(oclGpuContext, 1,
			(const char **) &oclSwScore, &program_length, &ciErrNum);
//	checkClError("Unable to build program.", ciErrNum);
	if (ciErrNum == CL_SUCCESS) {

		// build the program
		Log.Verbose("Build Options: %s", buildOptions.c_str());
		ciErrNum = clBuildProgram(cpProgram, 0, NULL, buildOptions.c_str(),
				NULL, NULL);
		if (ciErrNum != CL_SUCCESS)
				Log.Error("Build failed: %s", print_cl_errstring(ciErrNum));

		//checkClError("Unable to build program (clBuildProgram).", ciErrNum);
		//clUnloadCompiler();
		char cBuildLog[10240];

		clGetProgramBuildInfo(cpProgram, oclDevice, CL_PROGRAM_BUILD_OPTIONS,
				sizeof(cBuildLog), cBuildLog, NULL);
		Log.Verbose("Build options: %s", cBuildLog);

		cl_build_status status;
		clGetProgramBuildInfo(cpProgram, oclDevice, CL_PROGRAM_BUILD_STATUS,
				sizeof(status), &status, NULL);

		if (status != CL_BUILD_SUCCESS) {
			Log.Message("Build status: %s", print_cl_buildstatus(status));
			clGetProgramBuildInfo(cpProgram, oclDevice, CL_PROGRAM_BUILD_LOG,
					sizeof(cBuildLog), cBuildLog, NULL);

			Log.Message("Build log:");

			char * pBuildLog = strtok(cBuildLog, "\n");
			while (pBuildLog != NULL) {
				if (strlen(pBuildLog) > 1) {
					Log.Message("%s", pBuildLog);
				}
				pBuildLog = strtok(NULL, "\n");
			}
		}

		checkClError("Unable to build program end.", ciErrNum);
		return cpProgram;
	} else {
		Log.Error("Unable to load OpenCl kernel source. Error: %d", ciErrNum);
	}

	return 0;
}

const char* OclHost::print_cl_buildstatus(cl_build_status s) {
	switch (s) {
		case CL_BUILD_NONE: return "build not started";
		case CL_BUILD_ERROR: return "build failed";
		case CL_BUILD_SUCCESS: return "build succeeded";
		case CL_BUILD_IN_PROGRESS: return "build in progress";
		default: return "unknown";
	}
}

int OclHost::getThreadPerMulti() {
	if (isGPU()) {
#ifdef __APPLE__
		int id = 20;
#else
		cl_uint revision = getDeviceInfoInt(
		CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV);

		int id = 10 * revision
				+ getDeviceInfoLong(CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV);
#endif
		switch (id) {
		case 10:
			return 768;
		case 11:
			return 768;
		case 12:
			return 1024;
		case 13:
			return 1024;
		case 20:
			return 1536;
		case 21:
			return 1536;
		default:
			return 1536;
		}
	} else {
		return 256;
	}

}

//Maybe somebody could tell me how to use template when exporting a class from a dll. Probably not possible?
cl_uint OclHost::getDeviceInfoInt(cl_device_info info) {
	cl_uint value = 0;
	clGetDeviceInfo(oclDevice, info, sizeof(value), &value, 0);
	return value;
}

bool OclHost::isGPU() {
	return devType == CL_DEVICE_TYPE_GPU;
}

cl_ulong OclHost::getDeviceInfoLong(cl_device_info info) {
	cl_ulong value = 0;
	clGetDeviceInfo(oclDevice, info, sizeof(value), &value, 0);
	return value;
}

//cl_mem copytoDevice(size_t const size, cl_mem_flags flags, void * data) {
//	cl_int ciErrNum;
//	cl_mem gpu_var = clCreateBuffer(oclGpuContext, flags, size, data, &ciErrNum);
//	if (ciErrNum != CL_SUCCESS) {
//		Log.Error("Unable to create buffer on device. Error: %d", ciErrNum);
//	}
//	return gpu_var;
//}

void OclHost::writeToDevice(cl_mem buffer, cl_bool blocking_write,
		size_t offset, size_t size, const void * ptr) {

	cl_int ciErr = clEnqueueWriteBuffer(oclCommandQueue, buffer, blocking_write,
			offset, size, ptr, 0, NULL, NULL);
	clFlush(oclCommandQueue);
	//ciErr = clEnqueueWriteBuffer(oclCommandQueue[1], buffer, blocking_write, offset, size, ptr, 0, NULL, NULL);
	checkClError("Unable to write to Device.", ciErr);
}

void OclHost::readFromDevice(cl_mem buffer, cl_bool blocking_read,
		size_t offset, size_t size, void * ptr, size_t size_of) {
	cl_int ciErrNum = clEnqueueReadBuffer(oclCommandQueue, buffer,
			blocking_read, offset, size * size_of, ptr, 0, 0, 0);
	clFlush(oclCommandQueue);
	checkClError("Unable to read from device.", ciErrNum);
}

void OclHost::executeKernel(cl_kernel kernel, const size_t global_work_size,
		const size_t local_work_size) {
	cl_int ciErrNum;

	ciErrNum = clEnqueueNDRangeKernel(oclCommandQueue, kernel, 1, 0,
			&global_work_size, &local_work_size, 0, 0, 0);

	clFlush(oclCommandQueue);
	//ciErrNum = clEnqueueNDRangeKernel(oclCommandQueue[1], kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, 0);
	checkClError("Unable to execute kernel.", ciErrNum);
}

void OclHost::waitForDevice() {
	cl_int ciErr = clFinish(oclCommandQueue);
	//ciErr = clFinish(oclCommandQueue[1]);
	checkClError("clFinished failed.", ciErr);
}

char * OclHost::print_cl_errstring(cl_int err) {
	switch (err) {
	case CL_SUCCESS:
		return strdup("Success!");
	case CL_DEVICE_NOT_FOUND:
		return strdup("Device not found.");
	case CL_DEVICE_NOT_AVAILABLE:
		return strdup("Device not available");
	case CL_COMPILER_NOT_AVAILABLE:
		return strdup("Compiler not available");
	case CL_MEM_OBJECT_ALLOCATION_FAILURE:
		return strdup("Memory object allocation failure");
	case CL_OUT_OF_RESOURCES:
		return strdup("Out of resources");
	case CL_OUT_OF_HOST_MEMORY:
		return strdup("Out of host memory");
	case CL_PROFILING_INFO_NOT_AVAILABLE:
		return strdup("Profiling information not available");
	case CL_MEM_COPY_OVERLAP:
		return strdup("Memory copy overlap");
	case CL_IMAGE_FORMAT_MISMATCH:
		return strdup("Image format mismatch");
	case CL_IMAGE_FORMAT_NOT_SUPPORTED:
		return strdup("Image format not supported");
	case CL_BUILD_PROGRAM_FAILURE:
		return strdup("Program build failure");
	case CL_MAP_FAILURE:
		return strdup("Map failure");
	case CL_INVALID_VALUE:
		return strdup("Invalid value");
	case CL_INVALID_DEVICE_TYPE:
		return strdup("Invalid device type");
	case CL_INVALID_PLATFORM:
		return strdup("Invalid platform");
	case CL_INVALID_DEVICE:
		return strdup("Invalid device");
	case CL_INVALID_CONTEXT:
		return strdup("Invalid context");
	case CL_INVALID_QUEUE_PROPERTIES:
		return strdup("Invalid queue properties");
	case CL_INVALID_COMMAND_QUEUE:
		return strdup("Invalid command queue");
	case CL_INVALID_HOST_PTR:
		return strdup("Invalid host pointer");
	case CL_INVALID_MEM_OBJECT:
		return strdup("Invalid memory object");
	case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
		return strdup("Invalid image format descriptor");
	case CL_INVALID_IMAGE_SIZE:
		return strdup("Invalid image size");
	case CL_INVALID_SAMPLER:
		return strdup("Invalid sampler");
	case CL_INVALID_BINARY:
		return strdup("Invalid binary");
	case CL_INVALID_BUILD_OPTIONS:
		return strdup("Invalid build options");
	case CL_INVALID_PROGRAM:
		return strdup("Invalid program");
	case CL_INVALID_PROGRAM_EXECUTABLE:
		return strdup("Invalid program executable");
	case CL_INVALID_KERNEL_NAME:
		return strdup("Invalid kernel name");
	case CL_INVALID_KERNEL_DEFINITION:
		return strdup("Invalid kernel definition");
	case CL_INVALID_KERNEL:
		return strdup("Invalid kernel");
	case CL_INVALID_ARG_INDEX:
		return strdup("Invalid argument index");
	case CL_INVALID_ARG_VALUE:
		return strdup("Invalid argument value");
	case CL_INVALID_ARG_SIZE:
		return strdup("Invalid argument size");
	case CL_INVALID_KERNEL_ARGS:
		return strdup("Invalid kernel arguments");
	case CL_INVALID_WORK_DIMENSION:
		return strdup("Invalid work dimension");
	case CL_INVALID_WORK_GROUP_SIZE:
		return strdup("Invalid work group size");
	case CL_INVALID_WORK_ITEM_SIZE:
		return strdup("Invalid work item size");
	case CL_INVALID_GLOBAL_OFFSET:
		return strdup("Invalid global offset");
	case CL_INVALID_EVENT_WAIT_LIST:
		return strdup("Invalid event wait list");
	case CL_INVALID_EVENT:
		return strdup("Invalid event");
	case CL_INVALID_OPERATION:
		return strdup("Invalid operation");
	case CL_INVALID_GL_OBJECT:
		return strdup("Invalid OpenGL object");
	case CL_INVALID_BUFFER_SIZE:
		return strdup("Invalid buffer size");
	case CL_INVALID_MIP_LEVEL:
		return strdup("Invalid mip-map level");
	default:
		return strdup("Unknown");
	}
}
