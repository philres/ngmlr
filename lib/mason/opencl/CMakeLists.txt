cmake_minimum_required(VERSION 2.8)

ADD_EXECUTABLE(oclTool tools/oclTool.cpp)

IF(${CMAKE_MAJOR_VERSION} GREATER 2)
	ADD_CUSTOM_COMMAND(
	   OUTPUT ${PROJECT_BINARY_DIR}/oclDefines.h
	   COMMAND $<TARGET_FILE:oclTool> oclDefines ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclDefines.cl ${PROJECT_BINARY_DIR}/oclDefines.h
	   DEPENDS oclTool
	   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclDefines.cl
	)
	
	ADD_CUSTOM_COMMAND(
	   OUTPUT ${PROJECT_BINARY_DIR}/oclSwScore.h
	   COMMAND $<TARGET_FILE:oclTool> oclSwScore ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwScore.cl ${PROJECT_BINARY_DIR}/oclSwScore.h
	   DEPENDS oclTool
	   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwScore.cl
	)
	
	ADD_CUSTOM_COMMAND(
	   OUTPUT ${PROJECT_BINARY_DIR}/oclEndFreeScore.h
	   COMMAND $<TARGET_FILE:oclTool> oclEndFreeScore ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclEndFreeScore.cl ${PROJECT_BINARY_DIR}/oclEndFreeScore.h
	   DEPENDS oclTool
	   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclEndFreeScore.cl
	)
	
	ADD_CUSTOM_COMMAND(
	   OUTPUT ${PROJECT_BINARY_DIR}/oclSwAlignment.h
	   COMMAND $<TARGET_FILE:oclTool> oclSwAlignment ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwAlignment.cl ${PROJECT_BINARY_DIR}/oclSwAlignment.h
	   DEPENDS oclTool
	   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwAlignment.cl
	)
	
	ADD_CUSTOM_COMMAND(
	   OUTPUT ${PROJECT_BINARY_DIR}/oclSwCigar.h
	   COMMAND $<TARGET_FILE:oclTool> oclSwCigar ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwCigar.cl ${PROJECT_BINARY_DIR}/oclSwCigar.h
	   DEPENDS oclTool
	   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclSwCigar.cl
	)
else()
	GET_TARGET_PROPERTY(oclTool_exe oclTool LOCATION)
	ADD_CUSTOM_COMMAND(
	   OUTPUT ${PROJECT_BINARY_DIR}/oclDefines.h
	   COMMAND ${oclTool_exe} oclDefines ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclDefines.cl ${PROJECT_BINARY_DIR}/oclDefines.h
	   DEPENDS oclTool
	   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/lib/mason/opencl/opencl/oclDefines.cl
	)
	
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
endif()

include_directories (.)
include_directories (${PROJECT_BINARY_DIR})
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories( ${OPENCL_INCLUDE_DIRS} )

add_library(MASonOpenCl    STATIC ${PROJECT_BINARY_DIR}/oclDefines.h 
              ${PROJECT_BINARY_DIR}/oclSwScore.h
              ${PROJECT_BINARY_DIR}/oclEndFreeScore.h
              ${PROJECT_BINARY_DIR}/oclSwAlignment.h							
              ${PROJECT_BINARY_DIR}/oclSwCigar.h
              OclHost.cpp SWOcl.cpp SWOclAlignment.cpp SWOclCigar.cpp Timer.cpp)	

SET_TARGET_PROPERTIES(MASonOpenCl PROPERTIES COMPILE_FLAGS "-g -O3 -DNDEBUG")

IF(APPLE)
  target_link_libraries(MASonOpenCl "${OPENCL_LIBRARIES}")
ELSE()
  target_link_libraries(MASonOpenCl "-Wl,-Bdynamic ${OPENCL_LIBRARIES}")
ENDIF()

add_library(MASonOpenCl-debug    STATIC ${PROJECT_BINARY_DIR}/oclDefines.h 
              ${PROJECT_BINARY_DIR}/oclSwScore.h
              ${PROJECT_BINARY_DIR}/oclEndFreeScore.h
              ${PROJECT_BINARY_DIR}/oclSwAlignment.h							
              ${PROJECT_BINARY_DIR}/oclSwCigar.h
              OclHost.cpp SWOcl.cpp SWOclAlignment.cpp SWOclCigar.cpp Timer.cpp)	

set_target_properties(MASonOpenCl-debug PROPERTIES COMPILE_DEFINITIONS "DEBUGLOG")

SET_TARGET_PROPERTIES(MASonOpenCl-debug PROPERTIES COMPILE_FLAGS "-g3 -O0")

IF(APPLE)
  target_link_libraries(MASonOpenCl-debug "${OPENCL_LIBRARIES}")
ELSE()
  target_link_libraries(MASonOpenCl-debug "-Wl,-Bdynamic ${OPENCL_LIBRARIES}")
ENDIF()
