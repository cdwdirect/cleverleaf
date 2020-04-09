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
exec 3>&1 4>&2
function set_output_prefix() {
    exec 1>&3 2>&4
    exec > >(sed "s/^/$1-$2 [$3] ...........: /")
    exec 2> >(sed "s/^/$1-$2 [$3] (stderr) ..: /" >&2)
}

set_output_prefix $1 $2 "RAJA"
cd package-$1/$2/build/RAJA
make -j install
cd ../../../../

set_output_prefix $1 $2 "Umpire"
cd package-$1/$2/build/Umpire
make -j install || { echo "Failed building Umpire"; exit 1; }
cd ../../../../

set_output_prefix $1 $2 "SAMRAI"
cd package-$1/$2/build/SAMRAI
make -j install || { echo "Failed building SAMRAI"; exit 1; }
cd ../../../../
#
set_output_prefix $1 $2 "cleverleaf"
#./configure.$1 -m $2 -p package-$1 -c cleverleaf
cd package-$1/$2/build/cleverleaf
make -j install || { echo "Failed building cleverleaf"; exit 1; }
cd ../../../../
#
exec 1>&3 2>&4
echo "== CleverLeaf './go $1 $2' script is finished."
echo "== Binaries should be available in: ./package-$1/$2/install"
echo ""

