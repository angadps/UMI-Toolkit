#!/bin/sh
#$ -cwd
#$ -V
#$ -S /bin/bash

bamDir=/home/singhann/tools/bamtools-master
samDir=/home/singhann/tools/samtools-0.1.17

g++ -c -std=c++0x -O1 -g utils.cpp
g++ -c -std=c++0x -I $bamDir/include -I $samDir -O1 -g readUtils.cpp
g++ -c -std=c++0x -I $bamDir/include -I $samDir -O1 -g family.cpp
g++ -c -std=c++0x -I $bamDir/include -I $samDir -O1 -g main.cpp
g++ -c -std=c++0x -I $bamDir/include -I $samDir -O1 -g mutations.cpp

g++ -std=c++0x -L $bamDir/lib -L $samDir -O1 -g utils.o readUtils.o mutations.o family.o main.o -o main -lz -lbamtools  $samDir/libbam.a

#cp main ../

