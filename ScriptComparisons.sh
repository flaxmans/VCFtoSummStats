#!/bin/bash

git checkout master
make clean
make

./ConductSeveralTestRuns.sh 0

mv foo.out foo.master.7.27.out

git checkout buffering

make clean

make

echo "10k:"
./ConductSeveralTestRuns.sh 10000

mv foo.out foo.10kbuff.7.27.out

printf "\n\n100k:\n"
./ConductSeveralTestRuns.sh 100000

mv foo.out foo.100kbuff.7.27.out

printf "\n\n1M:\n"
./ConductSeveralTestRuns.sh 1000000

mv foo.out foo.1Mbuff.7.27.out

printf "\n\n10M:\n"
./ConductSeveralTestRuns.sh 10000000

mv foo.out foo.10Mbuff.7.27.out

printf "\n\n100M:\n"
./ConductSeveralTestRuns.sh 100000000

mv foo.out foo.100Mbuff.7.27.out

fnames="foo.10kbuff.7.27.out foo.100kbuff.7.27.out foo.1Mbuff.7.27.out foo.10Mbuff.7.27.out foo.100Mbuff.7.27.out"

printf "\n\nOutput files are:\n\t$fnames\n\n"
