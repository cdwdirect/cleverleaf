###############################################################################
# Copyright (c) 2016-19, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-689114
#
# All rights reserved.
#
# This file is part of RAJA.
#
# For details about use and distribution, please read RAJA/LICENSE.
#
###############################################################################

if (ENABLE_OPENMP)
  if(OPENMP_FOUND)
    list(APPEND RAJA_EXTRA_NVCC_FLAGS -Xcompiler ${OpenMP_CXX_FLAGS})
    message(STATUS "OpenMP Enabled")
  else()
    message(WARNING "OpenMP NOT FOUND")
    set(ENABLE_OPENMP Off)
  endif()
endif()

if (ENABLE_TBB)
  find_package(TBB)
  if(TBB_FOUND)
    blt_register_library(
      NAME tbb
      INCLUDES ${TBB_INCLUDE_DIRS}
      LIBRARIES ${TBB_LIBRARIES})
    message(STATUS "TBB Enabled")
  else()
    message(WARNING "TBB NOT FOUND")
    set(ENABLE_TBB Off)
  endif()
endif ()

if (ENABLE_CHAI)
  message(STATUS "CHAI enabled")
  find_package(umpire)
  find_package(chai)
  include_directories(${CHAI_INCLUDE_DIRS})
endif()

if (ENABLE_APOLLO)
  find_package(APOLLO REQUIRED)
  if(APOLLO_FOUND)
    blt_register_library(
      NAME apollo
      INCLUDES ${APOLLO_INCLUDE_DIRS}
      LIBRARIES ${APOLLO_LIBRARY})
    message(STATUS "APOLLO enabled")
    #message(STATUS "---- APOLLO:    APOLLO_INCLUDE_DIRS = ${APOLLO_INCLUDE_DIRS}")
    #message(STATUS "---- APOLLO:    APOLLO_LIB_DIRS     = ${APOLLO_LIB_DIRS}")
    #message(STATUS "---- APOLLO:    APOLLO_LIBRARY      = ${APOLLO_LIBRARY}")
    find_package(callpath REQUIRED)
    if(callpath_CONFIG_LOADED)
      message(STATUS "CALLPATH enabled")
      #message(STATUS "---- CALLPATH:    callpath_INCLUDE_DIR = ${callpath_INCLUDE_DIR}")
      #message(STATUS "---- CALLPATH:    callpath_LIB_DIR     = ${callpath_LIB_DIR}")
      blt_register_library(
        NAME callpath
        INCLUDES ${callpath_INCLUDE_DIR}
        LIBRARIES ${callpath_LIB_DIR}/libcallpath.so)
    else()
      message(STATUS "CALLPATH NOT FOUND")
      #message(STATUS "---- CALLPATH:  ** ERROR ** Could not locate the callpath library.")
    endif(callpath_CONFIG_LOADED)
  else()
    message(WARNING "APOLLO NOT FOUND")
  endif(APOLLO_FOUND)
endif()


