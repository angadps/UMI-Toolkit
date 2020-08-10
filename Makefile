#!/bin/sh
#$ -cwd
#$ -V
#$ -S /bin/bash

bamIDir=/home/singhann/tools/bamtools-master/include
bamLDir=/home/singhann/tools/bamtools-master/lib
samDir=/home/singhann/tools/samtools-0.1.17

CC=g++
CFLAGS= -std=c++0x
OFLAGS= -O1 -g

DEPS= utils.h readUtils.h mutations.h family.h defines.h
OBJ = utils.o readUtils.o mutations.o family.o build-consensus.o
LIBS= -lz -lbamtools

%.o: %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) -I $(bamIDir) -I $(samDir) $(OFLAGS) -o $@ $<

build-consensus: $(OBJ)
	$(CC) $(CFLAGS) -L $(bamLDir) -L $(samDir) $(OFLAGS) -o $@ $^ $(LIBS) $(samDir)/libbam.a

.PHONY: clean

clean:
	rm -f $(OBJ)

