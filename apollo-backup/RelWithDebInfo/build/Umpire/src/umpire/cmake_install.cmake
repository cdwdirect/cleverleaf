# Install script for directory: /g/g17/wood67/src/cleverleaf/Umpire/src/umpire

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/lib/libumpire.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/umpire" TYPE FILE FILES
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/Allocator.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/Replay.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/ResourceManager.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/ResourceManager.inl"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/TypedAllocator.hpp"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/TypedAllocator.inl"
    "/g/g17/wood67/src/cleverleaf/Umpire/src/umpire/Umpire.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/install/Umpire/include")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/install/Umpire" TYPE DIRECTORY FILES "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/include")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/src/umpire/tpl/cmake_install.cmake")
  include("/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/src/umpire/resource/cmake_install.cmake")
  include("/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/src/umpire/alloc/cmake_install.cmake")
  include("/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/src/umpire/op/cmake_install.cmake")
  include("/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/src/umpire/util/cmake_install.cmake")
  include("/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/Umpire/src/umpire/strategy/cmake_install.cmake")

endif()

