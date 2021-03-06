.. ###############################################################################
.. # Copyright (c) 2017, Lawrence Livermore National Security, LLC.
.. #
.. # Produced at the Lawrence Livermore National Laboratory
.. #
.. # LLNL-CODE-725085
.. #
.. # All rights reserved.
.. #
.. # This file is part of BLT.
.. #
.. # For additional details, please also read BLT/LICENSE.
.. #
.. # Redistribution and use in source and binary forms, with or without
.. # modification, are permitted provided that the following conditions are met:
.. #
.. # * Redistributions of source code must retain the above copyright notice,
.. #   this list of conditions and the disclaimer below.
.. #
.. # * Redistributions in binary form must reproduce the above copyright notice,
.. #   this list of conditions and the disclaimer (as noted below) in the
.. #   documentation and/or other materials provided with the distribution.
.. #
.. # * Neither the name of the LLNS/LLNL nor the names of its contributors may
.. #   be used to endorse or promote products derived from this software without
.. #   specific prior written permission.
.. #
.. # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
.. # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
.. # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
.. # ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
.. # LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
.. # DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
.. # DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
.. # OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
.. # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
.. # STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
.. # IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
.. # POSSIBILITY OF SUCH DAMAGE.
.. #
.. ###############################################################################

.. _Recommendations:

CMake Recommendations 
====================== 

This section includes several recommendations for how to wield CMake. Some of them are embodied in BLT, others are broader suggestions for CMake bliss.


.. rubric:: Disable in-source builds 

*BLT Enforces This*


In-source builds clutter source code with temporary build files and prevent other out-of-source builds from being created. Disabling in-source builds avoids clutter and accidental checkins of temporary build files.

.. rubric:: Avoid using globs to identify source files

Globs are evaluated at CMake configure time - not build time. This means CMake will not detect new source files when they are added to the file system unless there are other changes that trigger CMake to reconfigure. 

The CMake documentation also warns against this:
https://cmake.org/cmake/help/v3.10/command/file.html?highlight=glob#file


.. rubric::  Use arguments instead of options in CMake Macros and Functions

``CMAKE_PARSE_ARGUMENTS`` allows Macros or Functions to support options. Options are enabled by passing them by name when calling a Macro or Function. Because of this, wrapping an existing Macro or Function in a way that passes through options requires if tests and multiple copies of the call. For example:

.. code-block:: cmake

  if(OPTION)
      my_function(arg1 arg2 arg3 OPTION)
  else()
      my_function(arg1 arg2 arg3)
  endif()

Adding more options compounds the logic to achieve these type of calls.

To simplify calling logic, we recommend using an argument instead of an option.

.. code-block:: cmake

  if(OPTION)
      set(arg4_value ON)
  endif()
  
  my_function(arg1 arg2 arg3 ${arg4_value})


.. rubric::  Prefer explicit paths to locate third-party dependencies

Require passing explicit paths (ex: ``ZZZ_DIR``) for third-party dependency locations. This avoids surprises with incompatible installs sprinkled in various system locations. If you are using off-the-shelf *FindZZZ* logic, also consider adding CMake checks to verify that *FindZZZ* logic actually found the dependencies at the location specified.

.. rubric:: Emit a configure error if an explicitly identified third-party dependency is not found or an incorrect version is found.

If an explicit path to a dependency is given (ex: ``ZZZ_DIR``) it should be valid or result in a CMake configure error.

In contrast, if you only issue a warning and automatically disable a feature when a third-party dependency is bad, the warning often goes unnoticed and may not be caught until folks using your software are surprised. Emitting a configure error stops CMake and draws attention to the fact that something is wrong.  Optional dependencies are still supported by including them only if an explicit path to the dependency is given (ex: ``ZZZ_DIR``).



.. rubric::  Add headers as source files to targets

*BLT Macros Support This*

This ensures headers are tracked as dependencies and are included in the projects created by CMake's IDE generators, like Xcode or Eclipse. 


.. rubric::  Always support `make install`

This allows CMake to do the right thing based on ``CMAKE_INSTALL_PREFIX``, and also helps support CPack create release packages. This is especially important for libraries. In addition to targets, header files require an explicit install command.

Here is an example that installs a target and its headers:

.. code-block:: cmake

  ##################################
  # Install Targets for example lib
  ##################################
  install(FILES ${example_headers} DESTINATION include)
  install(TARGETS example
    EXPORT example
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
