cmake_minimum_required(VERSION 2.8)

ADD_EXECUTABLE(oclTool tools/oclTool.cpp)
GET_TARGET_PROPERTY(oclTool_exe oclTool LOCATION)

	#STATIC
    #set(CMAKE_EXE_LINKER_FLAGS -static) 
    #set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 
    #set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic 
    #set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS) 
    #set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC 
    #set(CMAKE_SHARED_LIBRARY_CXX_FLAGS) 
    #set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic 
    #set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS) 
    # Maybe this works as well, haven't tried yet. 
	#set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE) 	

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#SET (CMAKE_CXX_FLAGS " -g -O3 -Wall -fopenmp")
#SET (CMAKE_CXX_FLAGS_DEBUG "-D _DEBUG -g3 -O0 -Wall -fopenmp")

ADD_CUSTOM_COMMAND(
   OUTPUT ${PROJECT_BINARY_DIR}/oclDefines.h
   COMMAND ${oclTool_exe} oclDefines ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclDefines.cl ${PROJECT_BINARY_DIR}/oclDefines.h
   DEPENDS oclTool
   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclDefines.cl
)

message(${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwScore.cl)
ADD_CUSTOM_COMMAND(
   OUTPUT ${PROJECT_BINARY_DIR}/oclSwScore.h
   COMMAND ${oclTool_exe} oclSwScore ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwScore.cl ${PROJECT_BINARY_DIR}/oclSwScore.h
   DEPENDS oclTool
   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwScore.cl
)

ADD_CUSTOM_COMMAND(
   OUTPUT ${PROJECT_BINARY_DIR}/oclEndFreeScore.h
   COMMAND ${oclTool_exe} oclEndFreeScore ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclEndFreeScore.cl ${PROJECT_BINARY_DIR}/oclEndFreeScore.h
   DEPENDS oclTool
   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclEndFreeScore.cl
)

ADD_CUSTOM_COMMAND(
   OUTPUT ${PROJECT_BINARY_DIR}/oclSwAlignment.h
   COMMAND ${oclTool_exe} oclSwAlignment ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwAlignment.cl ${PROJECT_BINARY_DIR}/oclSwAlignment.h
   DEPENDS oclTool
   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwAlignment.cl
)

ADD_CUSTOM_COMMAND(
   OUTPUT ${PROJECT_BINARY_DIR}/oclSwCigar.h
   COMMAND ${oclTool_exe} oclSwCigar ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwCigar.cl ${PROJECT_BINARY_DIR}/oclSwCigar.h
   DEPENDS oclTool
   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwCigar.cl
)

include_directories (.)
include_directories (${PROJECT_BINARY_DIR})
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories( ${OPENCL_INCLUDE_DIRS} )
add_definitions(-fopenmp)


#add_library(MASonOpenClCpu SHARED ${PROJECT_BINARY_DIR}/oclDefines.h 
#							${PROJECT_BINARY_DIR}/oclSwScore.h
#							${PROJECT_BINARY_DIR}/oclEndFreeScore.h
#							${PROJECT_BINARY_DIR}/oclSwAlignment.h							
#							${PROJECT_BINARY_DIR}/oclSwCigar.h
#							SWOcl_export.cpp OclHost.cpp SWOcl.cpp SWOclAlignment.cpp SWOclCigar.cpp Timer.cpp)							

add_library(MASonOpenCl    STATIC ${PROJECT_BINARY_DIR}/oclDefines.h 
							${PROJECT_BINARY_DIR}/oclSwScore.h
							${PROJECT_BINARY_DIR}/oclEndFreeScore.h
							${PROJECT_BINARY_DIR}/oclSwAlignment.h							
							${PROJECT_BINARY_DIR}/oclSwCigar.h
							OclHost.cpp SWOcl.cpp SWOclAlignment.cpp SWOclCigar.cpp Timer.cpp)	

#add_library(MASonOpenClGpu SHARED ${PROJECT_BINARY_DIR}/oclDefines.h 
#							${PROJECT_BINARY_DIR}/oclSwScore.h
#							${PROJECT_BINARY_DIR}/oclEndFreeScore.h
#							${PROJECT_BINARY_DIR}/oclSwAlignment.h
#							${PROJECT_BINARY_DIR}/oclSwScore.h
#							${PROJECT_BINARY_DIR}/oclSwCigar.h
#							SWOcl_export.cpp OclHost.cpp  SWOcl.cpp SWOclAlignment.cpp SWOclCigar.cpp Timer.cpp)
							
#add_library(MASonOpenClGpu STATIC ${PROJECT_BINARY_DIR}/oclDefines.h 
#							${PROJECT_BINARY_DIR}/oclSwScore.h
#							${PROJECT_BINARY_DIR}/oclEndFreeScore.h
#							${PROJECT_BINARY_DIR}/oclSwAlignment.h
#							${PROJECT_BINARY_DIR}/oclSwScore.h
#							${PROJECT_BINARY_DIR}/oclSwCigar.h
#							SWOcl.cpp SWOclCigar.cpp SWOclAlignment.cpp Timer.cpp OclHost.cpp)							
							

#SET_TARGET_PROPERTIES(MASonOpenClCpu PROPERTIES
#                      DEFINE_SYMBOL CPU)

#target_link_libraries(MASonOpenClCpuS ${OPENCL_LIBRARIES})
#target_link_libraries(MASonOpenClGpuS ${OPENCL_LIBRARIES})

target_link_libraries(MASonOpenCl ${OPENCL_LIBRARIES})
#target_link_libraries(MASonOpenClGpu ${OPENCL_LIBRARIES})


ADD_CUSTOM_TARGET(opencl  
  DEPENDS
    MASonOpenClCpu MASonOpenCl 
)




