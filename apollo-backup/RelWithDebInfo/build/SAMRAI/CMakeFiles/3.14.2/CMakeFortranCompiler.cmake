set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-4.9.3/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "4.9.3")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")


set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/tce/packages/gcc/gcc-4.9.3/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3/finclude;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-dill-2.4-n44c2fod2h4gdfg7xdb5yhvldoqmnm2h/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-atl-2.2-sqsi7lqmh6rzj6xtxyk7izjhbbubf4vp/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libffs-1.5-6byc6slqd2k5xj7jw4sncuorndma2ux5/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libevpath-4.4.0-6eeomgq4j4hbkpe3pqncubqzmb7f5gya/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libnl-3.3.0-2cscqv3lmbamifzj72jm44j2mvb5j5ll/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-numpy-1.16.2-ckmegipscgcqnq2nrap3oxjzhxaoiqgd/lib/python2.7/site-packages/numpy/core/include;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.16-zsgpqf3kebcmbbrfxbtr2g5cf4pz2ymv/include;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3/include;/usr/local/include;/usr/tce/packages/gcc/gcc-4.9.3/include;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/x86_64-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/lib64;/lib64;/usr/lib64;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-dill-2.4-n44c2fod2h4gdfg7xdb5yhvldoqmnm2h/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/gtkorvo-atl-2.2-sqsi7lqmh6rzj6xtxyk7izjhbbubf4vp/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libffs-1.5-6byc6slqd2k5xj7jw4sncuorndma2ux5/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libevpath-4.4.0-6eeomgq4j4hbkpe3pqncubqzmb7f5gya/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/libnl-3.3.0-2cscqv3lmbamifzj72jm44j2mvb5j5ll/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-cffi-1.12.2-g7wwzrumfa54u6pphneeclelhotke3ga/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-numpy-1.16.2-ckmegipscgcqnq2nrap3oxjzhxaoiqgd/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-pandas-0.24.1-o2sxqnh3xwbc6n77h7tszg46q6ouhftj/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-scikit-learn-0.20.2-34kpy2ueayabksqehbcaeks3muij2xij/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/py-scipy-1.2.1-7qvyjedp4mphsyfjqmemgdeqfv5ls4vz/lib;/g/g17/wood67/src/spack/opt/spack/linux-rhel7-x86_64/gcc-4.9.3/python-2.7.16-zsgpqf3kebcmbbrfxbtr2g5cf4pz2ymv/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
