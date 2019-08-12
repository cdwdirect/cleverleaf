set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "4.9.3")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "90")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc-ranlib")
set(CMAKE_LINKER "/usr/tce/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-dill-2.4-n44c2fod2h4gdfg7xdb5yhvldoqmnm2h/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-atl-2.2-sqsi7lqmh6rzj6xtxyk7izjhbbubf4vp/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libffs-1.5-6byc6slqd2k5xj7jw4sncuorndma2ux5/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libevpath-4.4.0-6eeomgq4j4hbkpe3pqncubqzmb7f5gya/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libnl-3.3.0-2cscqv3lmbamifzj72jm44j2mvb5j5ll/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-numpy-1.16.2-ckmegipscgcqnq2nrap3oxjzhxaoiqgd/lib/python2.7/site-packages/numpy/core/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.16-zsgpqf3kebcmbbrfxbtr2g5cf4pz2ymv/include;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3/include;/usr/local/include;/usr/tce/packages/gcc/gcc-4.9.3/include;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3/include-fixed;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/lib64;/lib64;/usr/lib64;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-dill-2.4-n44c2fod2h4gdfg7xdb5yhvldoqmnm2h/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-atl-2.2-sqsi7lqmh6rzj6xtxyk7izjhbbubf4vp/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libffs-1.5-6byc6slqd2k5xj7jw4sncuorndma2ux5/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libevpath-4.4.0-6eeomgq4j4hbkpe3pqncubqzmb7f5gya/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libnl-3.3.0-2cscqv3lmbamifzj72jm44j2mvb5j5ll/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-cffi-1.12.2-g7wwzrumfa54u6pphneeclelhotke3ga/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-numpy-1.16.2-ckmegipscgcqnq2nrap3oxjzhxaoiqgd/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-pandas-0.24.1-o2sxqnh3xwbc6n77h7tszg46q6ouhftj/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-scikit-learn-0.20.2-34kpy2ueayabksqehbcaeks3muij2xij/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-scipy-1.2.1-7qvyjedp4mphsyfjqmemgdeqfv5ls4vz/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.16-zsgpqf3kebcmbbrfxbtr2g5cf4pz2ymv/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
