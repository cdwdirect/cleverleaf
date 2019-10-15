#!/bin/bash
if [ "x$1" == "x" ] ; then
    echo "USAGE: ./go <apollo | normal> <Release | RelWithDebInfo>"
    kill -INT $$
fi
if [ "x$2" == "x" ] ; then
    echo "USAGE: ./go <apollo | normal> <Release | RelWithDebInfo>"
    kill -INT $$
fi
#
cd package-$1/$2/build/RAJA
make -j install
cd ../../../../
#
cd package-$1/$2/build/Umpire
make -j install
cd ../../../../
#
cd package-$1/$2/build/SAMRAI
make -j install
cd ../../../../
#
cd package-$1/$2/build/cleverleaf
make -j install
cd ../../../../
#
echo "== CleverLeaf './go $1 $2' script is finished."
echo "== Binaries should be available in: ./package-$1/$2/install"
echo ""

