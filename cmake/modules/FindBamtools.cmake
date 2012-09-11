# - Try to find Bamtools
# Once done this will define
#  BAMTOOLS_FOUND - System has LibXml2
#  BAMTOOLS_INCLUDE_DIRS - The LibXml2 include directories
#  BAMTOOLS_LIBRARIES - The libraries needed to use LibXml2

#find_path(BamTools_INCLUDE_DIR api/api_global.h )
#find_library(BamTools_LIBRARY NAMES bamtools HINTS ${BamTools_PREFIX}/lib/bamtools )

find_path(BAMTOOLS_INCLUDE_DIR api/api_global.h QUIET)
find_library(BAMTOOLS_LIBRARY NAMES libbamtools.a)

set(BAMTOOLS_LIBRARIES ${BAMTOOLS_LIBRARY} )
set(BAMTOOLS_INCLUDE_DIRS ${BAMTOOLS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BamTools  DEFAULT_MSG
                                  BAMTOOLS_LIBRARY BAMTOOLS_INCLUDE_DIR)

mark_as_advanced(BAMTOOLS_INCLUDE_DIR BAMTOOLS_LIBRARY )