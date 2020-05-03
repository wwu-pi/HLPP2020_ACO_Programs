#!/bin/bash

# remove files and create folder
rm -rf -- build && \
mkdir build && \

# run cmake
cd build && \
cmake -G "Unix Makefiles" -D CMAKE_BUILD_TYPE=dev -D CMAKE_CXX_COMPILER=g++ ../ && \

make aco_iroulette_0 && \
echo "Runs; Iterations;Problem;ants;BestSequence;RouteDistance;Time;"
for CITY in 2 4
do
	for ANT in 1024 2048 4096 8192
	do
		for RUN in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
		do
			bin/aco_iroulette_0 $RUN 15 $CITY $ANT
		done
	done
done
