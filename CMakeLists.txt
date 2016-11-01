# ===========================================================================
#                                 SimBA-hap
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (SimBA-hap CXX)
set (SEQAN_APP_VERSION "0.0.1")

# ----------------------------------------------------------------------------
# Project Configuration
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake release mode." FORCE)
endif ()

# Set output binary directory.
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# ----------------------------------------------------------------------------
# Load SeqAn Build System
# ----------------------------------------------------------------------------

message (STATUS "Configuring SeqAn")
set (CMAKE_INCLUDE_PATH
    ${CMAKE_SOURCE_DIR}/../seqan/include
    ${CMAKE_INCLUDE_PATH})
set (CMAKE_MODULE_PATH
#    ${CMAKE_SOURCE_DIR}/../seqan/util/cmake/FindTBB
    ${CMAKE_SOURCE_DIR}/../seqan/util/cmake
    ${CMAKE_MODULE_PATH})

include (SeqAnContribs)
include (SeqAnBuildSystem)
set (SEQAN_ROOT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../seqan/" CACHE INTERNAL "Root source directory." FORCE)
seqan_build_system_init ()
seqan_get_repository_info ()
seqan_setup_library ()
include (package)

# ----------------------------------------------------------------------------
# Include Sources Directory
# ----------------------------------------------------------------------------

add_subdirectory (src)

# ----------------------------------------------------------------------------
# CPack Configuration
# ----------------------------------------------------------------------------

set (CPACK_PACKAGE_NAME "SimBA-hap")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Phasing")
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <esiragu@us.ibm.com>")
set (CPACK_PACKAGE_VENDOR "Enrico Siragusa <esiragu@us.ibm.com>")

set (CPACK_PACKAGE_VERSION "0.0.1")
set (CPACK_PACKAGE_VERSION_MAJOR "0")
set (CPACK_PACKAGE_VERSION_MINOR "0")
set (CPACK_PACKAGE_VERSION_PATCH "1")

set (CPACK_GENERATOR "TGZ")

include (CPack)
