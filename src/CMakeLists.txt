# ===========================================================================
#                                 SimBA-hap
# ===========================================================================

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package (ZLIB)
find_package (SeqAn REQUIRED)
#find_package (Boost 1.51.0 REQUIRED)

if (NOT CXX11_FOUND)
    message (STATUS "  C++11 not found: built aborted")
    return ()
endif (NOT CXX11_FOUND)

#if (NOT SeqAn_FOUND)
#  message (STATUS "  SeqAn library not found: build aborted")
#  return ()
#endif (NOT SeqAn_FOUND)

#if (NOT Boost_FOUND)
#  message (STATUS "  Boost library not found: build aborted")
#  return ()
#endif (NOT Boost_FOUND)

if (NOT ZLIB_FOUND)
  message (STATUS "  Zlib library not found: gzip support disabled")
endif (NOT ZLIB_FOUND)


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories for SeqAn and other dependencies.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${LEMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by find_package.
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Apps
# ----------------------------------------------------------------------------

# Add executables.
add_executable (sim sim.cpp)
target_link_libraries (sim ${SEQAN_LIBRARIES} ${LEMON_LIBRARIES})

# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install in ${PREFIX}/bin directory
#install (TARGETS sim DESTINATION bin)
