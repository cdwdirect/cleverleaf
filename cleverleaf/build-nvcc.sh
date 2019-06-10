#!/bin/bash

. /usr/local/tools/dotkit/init.sh

use cmake

. /g/g92/ukbeck/Projects/spack/share/spack/setup-env.sh

module load cudatoolkit/7.5

spack use SAMRAI%gcc
spack use openmpi%gcc
spack use hdf5%gcc

export LIBRARY_PATH=""
export CMAKE_PREFIX_PATH=/g/g92/ukbeck/Projects/cudapatchdata/chaos_5_x86_64_ib:$CMAKE_PREFIX_PATH
export RAJA_DIR=/g/g92/ukbeck

rm -rf ${SYS_TYPE}_nvcc
mkdir ${SYS_TYPE}_nvcc

cd ${SYS_TYPE}_nvcc
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCLEVERLEAF_USE_GPU=On \
  ..

make -j 16 VERBOSE=1
