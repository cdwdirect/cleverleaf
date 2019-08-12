#!/bin/bash
if [ "x$1" == "x" ] ; then
    echo "USAGE: ./go <apollo | normal>"
    kill -INT $$
fi
#
cd package-$1/RelWithDebInfo/build/RAJA
make -j install
cd ../../../../
#
cd package-$1/RelWithDebInfo/build/Umpire
make -j install
cd ../../../../
#
cd package-$1/RelWithDebInfo/build/SAMRAI
make -j install
cd ../../../../
#
cd package-$1/RelWithDebInfo/build/cleverleaf
make -j install
cd ../../../../
#
echo "== CleverLeaf './go $1' script is finished."
echo "== Binaries should be available in: ./package-$1/RelWithDebInfo/install"
echo ""

