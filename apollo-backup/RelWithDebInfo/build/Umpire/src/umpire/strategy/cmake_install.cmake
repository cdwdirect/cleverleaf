# Install script for directory: /g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/umpire/strategy" TYPE FILE FILES
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/AllocationAdvisor.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/AllocationStrategy.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/AllocationTracker.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/DynamicPool.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/DynamicPoolHeuristic.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/FixedPool.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/MixedPool.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/MonotonicAllocationStrategy.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/SlotPool.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/SizeLimiter.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/ThreadSafeAllocator.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/umpire/strategy/mixins" TYPE FILE FILES "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/strategy/mixins/Inspector.hpp")
endif()

