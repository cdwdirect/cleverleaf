#!/bin/bash
cd apollo-test/RelWithDebInfo/build/RAJA
make -j install
cd ../../../../
cd apollo-test/RelWithDebInfo/build/Umpire
make -j install
cd ../../../../
cd apollo-test/RelWithDebInfo/build/SAMRAI
make -j install
cd ../../../../
cd apollo-test/RelWithDebInfo/build/cleverleaf
make -j install
cd ../../../../

