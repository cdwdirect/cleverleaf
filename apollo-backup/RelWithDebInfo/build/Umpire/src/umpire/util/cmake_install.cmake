# Install script for directory: /g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/install/Umpire")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/umpire/util" TYPE FILE FILES
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/AllocationMap.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/AllocationRecord.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/allocation_statistics.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/detect_vendor.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/Exception.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/Logger.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/Macros.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/MemoryResourceTraits.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/util/Platform.hpp"
    )
endif()

