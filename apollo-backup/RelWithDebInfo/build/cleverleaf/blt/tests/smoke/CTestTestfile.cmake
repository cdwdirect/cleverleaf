# CMake generated Testfile for 
# Source directory: /g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke
# Build directory: /g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/cleverleaf/blt/tests/smoke
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(blt_gtest_smoke "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/cleverleaf/tests/blt_gtest_smoke")
set_tests_properties(blt_gtest_smoke PROPERTIES  _BACKTRACE_TRIPLES "/g/g17/wood67/src/cleverleaf/cleverleaf/blt/cmake/BLTMacros.cmake;709;add_test;/g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke/CMakeLists.txt;57;blt_add_test;/g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke/CMakeLists.txt;0;")
add_test(blt_openmp_smoke "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/cleverleaf/tests/blt_openmp_smoke")
set_tests_properties(blt_openmp_smoke PROPERTIES  _BACKTRACE_TRIPLES "/g/g17/wood67/src/cleverleaf/cleverleaf/blt/cmake/BLTMacros.cmake;709;add_test;/g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke/CMakeLists.txt;116;blt_add_test;/g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke/CMakeLists.txt;0;")
add_test(blt_mpi_smoke "/usr/bin/srun" "-n" "2" "/g/g17/wood67/src/cleverleaf/apollo-test/RelWithDebInfo/build/cleverleaf/tests/blt_mpi_smoke")
set_tests_properties(blt_mpi_smoke PROPERTIES  _BACKTRACE_TRIPLES "/g/g17/wood67/src/cleverleaf/cleverleaf/blt/cmake/BLTMacros.cmake;709;add_test;/g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke/CMakeLists.txt;131;blt_add_test;/g/g17/wood67/src/cleverleaf/cleverleaf/blt/tests/smoke/CMakeLists.txt;0;")
