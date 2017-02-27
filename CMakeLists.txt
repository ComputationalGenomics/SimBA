# ===========================================================================
#                                 SimBA-hap
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (SimBA-hap CXX)
set (SEQAN_APP_VERSION "1.0.0")

# ----------------------------------------------------------------------------
# Project Configuration
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake release mode." FORCE)
endif ()

option (STATIC_BUILDS "Build static binaries." OFF)

# Set output binary directory.
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Enable static linkage.
if (STATIC_BUILDS)
    message (STATUS "Building static binaries.")

    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    if (APPLE)
        # Static build not supported on apple, but at least we can include gcc libs
        if (COMPILER_GCC)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        endif (COMPILER_GCC)
    else (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

        # Make sure -rdynamic isn't added automatically.
        set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
        # Make sure -fPIC isn't added automatically.
        set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)

        # For unknown reasons finding .a only seems to work for libz and
        # libbzip2; cmake than proceeds to wrap these in
        # -Wl,-Bstatic -lz -lbz2 -Wl,-Bdynamic
        # the latter reactivates dynamic linking for the system libs
        # we override this behaviour here:
        set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    endif (APPLE)
endif ()

# Strip binaries when packaging.
if ((CMAKE_BUILD_TYPE STREQUAL "Release") AND
    (NOT APPLE) AND
    (COMPILER_CLANG OR COMPILER_GCC))
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")
endif ()

# ----------------------------------------------------------------------------
# Load SeqAn Build System
# ----------------------------------------------------------------------------

message (STATUS "Configuring SeqAn")
set (CMAKE_INCLUDE_PATH
    ${CMAKE_SOURCE_DIR}/deps/seqan/include
    ${CMAKE_INCLUDE_PATH})
set (CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/deps/seqan/util/cmake
    ${CMAKE_MODULE_PATH})
set (CMAKE_PREFIX_PATH
    ${CMAKE_SOURCE_DIR}/deps/seqan/util/cmake
    ${CMAKE_PREFIX_PATH})

# ----------------------------------------------------------------------------
# Prepare Lemon Library
# ----------------------------------------------------------------------------

message (STATUS "Configuring Lemon")

FIND_PATH(LEMON_SOURCE_ROOT_DIR CMakeLists.txt
  PATHS ${CMAKE_SOURCE_DIR}/deps/lemon
  NO_DEFAULT_PATH
  DOC "Location of LEMON source as a CMAKE subproject")

IF(EXISTS ${LEMON_SOURCE_ROOT_DIR})
  ADD_SUBDIRECTORY(${LEMON_SOURCE_ROOT_DIR} deps/lemon)
  SET(LEMON_INCLUDE_DIRS
    ${LEMON_SOURCE_ROOT_DIR}
    ${CMAKE_BINARY_DIR}/deps/lemon
  )
  SET(LEMON_LIBRARIES lemon)
  UNSET(LEMON_ROOT_DIR CACHE)
  UNSET(LEMON_DIR CACHE)
  UNSET(LEMON_INCLUDE_DIR CACHE)
  UNSET(LEMON_LIBRARY CACHE)
ELSE()
  FIND_PACKAGE(LEMON QUIET NO_MODULE)
  FIND_PACKAGE(LEMON REQUIRED)
ENDIF()

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
set (CPACK_PACKAGE_VENDOR "IBM Research")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "simba-hap ${CPACK_PACKAGE_VERSION}")

# Install non-binary files in "." directory.
install (FILES LICENSE README.md DESTINATION .)

if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(CPACK_GENERATOR "ZIP;DragNDrop")
elseif (CMAKE_SYSTEM_NAME MATCHES "Linux")
  set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB;RPM")
endif ()

set (CPACK_PACKAGE_VERSION "1.0.0")
set (CPACK_PACKAGE_VERSION_MAJOR "1")
set (CPACK_PACKAGE_VERSION_MINOR "0")
set (CPACK_PACKAGE_VERSION_PATCH "0")

set (CPACK_GENERATOR "TGZ")

include (CPack)
