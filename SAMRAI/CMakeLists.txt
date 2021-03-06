cmake_minimum_required(VERSION 3.1)

project(SAMRAI C CXX Fortran)

set(SAMRAI_VERSION_MAJOR 3)
set(SAMRAI_VERSION_MINOR 11)
set(SAMRAI_VERSION_PATCHLEVEL 1)
set(SAMRAI_VERSION
  "${SAMRAI_VERSION_MAJOR}.${SAMRAI_VERSION_MINOR}.${SAMRAI_VERSION_PATCHLEVEL}")


cmake_policy(SET CMP0057 NEW)

if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
  message(FATAL_ERROR "The BLT submodule is not present. \
If in the git repository, please run the following commands:\n \
git submodule init\n \
git submodule update")
endif ()

# BLT Options need to be set in the CACHE
set(ENABLE_FORTRAN On CACHE Bool "Enable Fortran")
set(ENABLE_OPENMP On CACHE Bool "Enable OpenMP")
set(ENABLE_CUDA On CACHE Bool "Enable CUDA")
set(ENABLE_MPI On CACHE Bool "Enable MPI.")

# Third party libraries
option(ENABLE_APOLLO "Enable Apollo" ON)
option(ENABLE_UMPIRE "Enable UMPIRE memory management" ON)
option(ENABLE_RAJA "Enable RAJA" ON)
option(ENABLE_CONDUIT "Enable CONDUIT" Off)
option(ENABLE_HDF5 "Enable HDF5." On)
option(ENABLE_HYPRE "Enable HYPRE" Off)
option(ENABLE_PETSC "Enable PETSc" Off)
option(ENABLE_PTSCOTCH "Enable Scotch" Off)
option(ENABLE_SILO "Enable Silo" Off)
option(ENABLE_SUNDIALS "Enable SUNDIALS" Off)

# SAMRAI options
option(ENABLE_CHECK_ASSERTIONS "Enable assertion checking." Off)
option(ENABLE_CHECK_DEV_ASSERTIONS "Enable SAMRAI developer assertion checking." Off)
option(ENABLE_CHECK_DIM_ASSERTIONS "Enable assertion checking for dimensions." Off)
option(ENABLE_BOX_COUNTING "Turns on box telemetry." Off)
option(ENABLE_DEPRECATED "Build with deprecated features." Off)
option(ENABLE_TIMERS "Enable SAMRAI timers." Off)
option(DEBUG_INITIALIZE_UNDEFINED "Initialize new memory to undefined values" Off)
set(MAXDIM 3 CACHE INT "Maximum allowed dimension")
option(ENABLE_NVTX_REGIONS "Enable NVTX regions with push/pop for NVVP." Off)
set(CUDA_ARCH "sm_35" CACHE STRING "Compute architecture to pass to CUDA builds")
set(CMAKE_CUDA_FLAGS "" CACHE STRING "")

include(GNUInstallDirs)

include(blt/SetupBLT.cmake)

# Set specific options for CUDA if enabled
if (ENABLE_RAJA AND ENABLE_CUDA)
  # RAJA requires some experimental features
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
endif ()

# TODO Ensure openmp flags are not enabled twice!
if (ENABLE_OPENMP AND ENABLE_CUDA)
  # CMAKE_CUDA_FLAGS does not pass OpenMP_CXX_FLAGS to the host compiler by default
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${OpenMP_CXX_FLAGS}")
endif ()

if (NOT ENABLE_UMPIRE)
  message(FATAL_ERROR "This code branch requires UMPIRE, for now...")
endif ()

if (ENABLE_CUDA)
  if (NOT ENABLE_RAJA)
    message(FATAL_ERROR "CUDA support requires RAJA")
  endif ()

  if (NOT ENABLE_UMPIRE)
    message(FATAL_ERROR "CUDA support requires UMPIRE")
  endif ()
endif ()

find_package(callpath REQUIRED)
if(callpath_CONFIG_LOADED)
    message(STATUS "---- CALLPATH:  Found!")
    message(STATUS "---- CALLPATH:    callpath_INCLUDE_DIR = ${callpath_INCLUDE_DIR}")
    message(STATUS "---- CALLPATH:    callpath_LIB_DIR     = ${callpath_LIB_DIR}")
    blt_register_library(
        NAME callpath
        INCLUDES ${callpath_INCLUDE_DIR}
        LIBRARIES ${callpath_LIB_DIR}/libcallpath.so)
else()
    message(STATUS "---- CALLPATH:  ** ERROR ** Could not locate the callpath library.")
endif()


if (ENABLE_OPENMP_TRACE)
    add_definitions(" -DRAJA_ENABLE_OPENMP_TRACE -DENABLE_OPENMP_TRACE ")
endif()

#if (ENABLE_APOLLO)
#  message(STATUS "Building SAMRAI w/APOLLO Support")
#  find_package(APOLLO REQUIRED)
#  if(APOLLO_FOUND)
#      add_definitions("-DENABLE_APOLLO")
#      blt_register_library(
#      NAME apollo
#      INCLUDES ${APOLLO_INCLUDE_DIRS}
#      LIBRARIES ${APOLLO_LIBRARY})
#      message(STATUS "---- APOLLO:    APOLLO_INCLUDE_DIRS = ${APOLLO_INCLUDE_DIRS}")
#      message(STATUS "---- APOLLO:    APOLLO_LIB_DIRS     = ${APOLLO_LIB_DIRS}")
#      message(STATUS "---- APOLLO:    APOLLO_LIBRARY      = ${APOLLO_LIBRARY}")
#  endif(APOLLO_FOUND)
#endif ()

include_directories(${CMAKE_BINARY_DIR}/include)

include(cmake/CMakeBasics.cmake)

set(SAMRAI_LIBRARIES SAMRAI_algs SAMRAI_appu SAMRAI_geom SAMRAI_hier
    SAMRAI_math SAMRAI_mesh SAMRAI_pdat SAMRAI_solv SAMRAI_tbox
    SAMRAI_xfer)

add_subdirectory(source)
add_subdirectory(tools)
add_subdirectory(docs)

export(TARGETS ${SAMRAI_LIBRARIES}
  FILE ${CMAKE_CURRENT_BINARY_DIR}/SAMRAITargets.cmake)

set(SAMRAI_MODULE "share/SAMRAI/cmake")
file(RELATIVE_PATH REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}" "${CMAKE_INSTALL_PREFIX}/include")

set(CONF_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/${REL_INCLUDE_DIR}")
configure_file(SAMRAIConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfig.cmake" @ONLY)


configure_file(SAMRAIConfigVersion.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfigVersion.cmake" @ONLY)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${SAMRAI_MODULE}
  )

install(EXPORT SAMRAITargets DESTINATION
  ${CMAKE_INSTALL_PREFIX}/${SAMRAI_MODULE})

install(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/include/SAMRAI
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
