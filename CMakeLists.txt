# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/vcf_sim/CMakeLists.txt
#
# CMakeLists.txt file for vcf_sim.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_vcf_sim CXX)
message (STATUS "Configuring apps/vcf_sim")

set (SEQAN_APP_VERSION "0.0.1")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
find_package (ZLIB)
find_package (SeqAn REQUIRED)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (STATUS "  Zlib not found: building vcf_sim without gzipped input/output.")
endif (NOT ZLIB_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (vcf_sim vcf_sim.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (vcf_sim ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (vcf_sim)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install vcf_sim in ${PREFIX}/bin directory
install (TARGETS vcf_sim
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/vcf_sim for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (vcf_sim)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:vcf_sim")
  set (CPACK_PACKAGE_NAME "vcf_sim")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "vcf_sim")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <esiragu@us.ibm.com>")
  set (CPACK_PACKAGE_VENDOR "Enrico Siragusa <esiragu@us.ibm.com>")

  seqan_configure_cpack_app (vcf_sim "vcf_sim")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:vcf_sim")

