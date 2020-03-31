cmake_minimum_required (VERSION 3.1)
project (CleverLeaf CXX)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
  message(FATAL_ERROR "The BLT submodule is not present. \
If in the git repository, please run the following commands:\n \
git submodule init\n \
git submodule update")
endif()

set(ENABLE_APOLLO ON CACHE Bool "Enable Apollo")
set(ENABLE_MPI ON CACHE BOOL "")
set(ENABLE_OPENMP ON CACHE Bool "Enable OpenMP")
set(ENABLE_CUDA ON CACHE Bool "Enable CUDA")

set (BLT_CXX_STD "c++14")
include(blt/SetupBLT.cmake)

if(ENABLE_CUDA)
  set(CMAKE_CUDA_FLAGS "--expt-extended-lambda --expt-relaxed-constexpr -DENABLE_UVM ${CMAKE_CUDA_FLAGS}")
endif()

if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        message(STATUS "---- OPENMP:  Found!")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else()
        message(STATUS "---- OPENMP:  Not found.")
    endif(OPENMP_FOUND)
endif(ENABLE_OPENMP)

if (ENABLE_APOLLO)
  message(STATUS "Building cleverleaf w/APOLLO Support")
  find_package(APOLLO REQUIRED)
  if(APOLLO_FOUND)
      add_definitions("-DENABLE_APOLLO")
    blt_register_library(
      NAME apollo
      INCLUDES ${APOLLO_INCLUDE_DIRS}
      LIBRARIES ${APOLLO_LIBRARY})
      message(STATUS "---- APOLLO:    APOLLO_INCLUDE_DIRS = ${APOLLO_INCLUDE_DIRS}")
      message(STATUS "---- APOLLO:    APOLLO_LIB_DIRS     = ${APOLLO_LIB_DIRS}")
      message(STATUS "---- APOLLO:    APOLLO_LIBRARY      = ${APOLLO_LIBRARY}")
  endif(APOLLO_FOUND)
endif ()

find_package(callpath REQUIRED)
if(callpath_CONFIG_LOADED)
  message(STATUS "CALLPATH enabled")
  message(STATUS "---- CALLPATH:    callpath_INCLUDE_DIR = ${callpath_INCLUDE_DIR}")
  message(STATUS "---- CALLPATH:    callpath_LIB_DIR     = ${callpath_LIB_DIR}")
  blt_register_library(
    NAME callpath
    INCLUDES ${callpath_INCLUDE_DIR}
    LIBRARIES ${callpath_LIB_DIR}/libcallpath.so)
else()
  message(STATUS "CALLPATH NOT FOUND")
endif(callpath_CONFIG_LOADED)

# ----- deps -----
find_package(RAJA REQUIRED)
find_package(umpire REQUIRED)
find_package(HDF5 REQUIRED)
find_package(SAMRAI REQUIRED)

# ----- system and build definitions -----
# Get the host name
cmake_host_system_information(RESULT SYSINFO QUERY HOSTNAME)
list(GET SYSINFO 0, HOST_NAME)

# Get the current working branch
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  add_definitions(-DVERSION="${GIT_BRANCH}" -DHOST_NAME="${HOST_NAME}")
endif ()

add_subdirectory (src)
