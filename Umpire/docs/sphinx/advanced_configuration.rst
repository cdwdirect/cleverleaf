.. _advanced_configuration:

======================
Advanced Configuration
======================

In addition to the normal options provided by CMake, Umpire uses some additional
configuration arguments to control optional features and behavior. Each
argument is a boolean option, and  can be turned on or off:

.. code-block:: bash

    -DENABLE_CUDA=Off

Here is a summary of the configuration options, their default value, and meaning:

      ===========================  ======== ===============================================================================
      Variable                     Default  Meaning
      ===========================  ======== ===============================================================================
      ``ENABLE_CUDA``              On       Enable CUDA support
      ``ENABLE_NUMA``              Off      Enable NUMA support
      ``ENABLE_HCC``               Off      Enable HCC support
      ``ENABLE_STATISTICS``        Off      Enable collection of memory statistics
      ``ENABLE_TESTING``           On       Build test executables
      ``ENABLE_BENCHMARKS``        On       Build benchmark programs
      ``ENABLE_LOGGING``           On       Enable Logging within Umpire
      ``ENABLE_SLIC``              Off      Enable SLIC logging
      ``ENABLE_ASSERTS``           On       Enable UMPIRE_ASSERT() within Umpire
      ``ENABLE_TOOLS``             On       Enable tools like replay
      ``ENABLE_DOCS``              Off      Build documentation (requires Sphinx and/or Doxygen)
      ``ENABLE_C``             Off      Build the C API
      ``ENABLE_FORTRAN``       Off      Build the Fortran API
      ===========================  ======== ===============================================================================

These arguments are explained in more detail below:

* ``ENABLE_CUDA``
  This option enables support for NVIDIA GPUs. If Umpire is built without CUDA
  or HCC support, then only the ``HOST`` allocator is available for use.

* ``ENABLE_NUMA``
  This option enables support for NUMA. The
  :class:`umpire::strategy::NumaPolicy` is available when built with this
  option, which may be used to locate the allocation to a specific node.

* ``ENABLE_HCC``
  This option enables support for AMD GPUs using the ROCm stack and HCC
  programming model. If Umpire is built without CUDA or HCC support, then only
  the ``HOST`` allocator is available for use.

* ``ENABLE_STATISTICS``
  This option enables collection of memory statistics. If Umpire is built with
  this option, the Conduit library will also be built.

* ``ENABLE_TESTING``
  This option controls whether or not test executables will be built.

* ``ENABLE_BENCHMARKS``
  This option will build the benchmark programs used to test performance.

* ``ENABLE_LOGGING``
  This option enables usage of Logging services for Umpire

* ``ENABLE_SLIC``
  This option enables usage of logging services provided by SLIC.

* ``ENABLE_ASSERTS``
  Enable assert() within Umpire

* ``ENABLE_TOOLS``
  Enable development tools for Umpire (replay, etc.)

* ``ENABLE_DOCS``
  Build user documentation (with Sphinx) and code documentation (with Doxygen)

* ``ENABLE_C``
  Build the C API, this allows accessing Umpire Allocators and the
  ResourceManager through a C interface.

* ``ENABLE_FORTRAN``
  Build the Fortran API.
