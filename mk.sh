#!/bin/bash

echo "-------------------------------------------"
echo "|                PntWrks                  |"
echo "-------------------------------------------"
mkdir sim/out
mkdir mod
rm -r sim/out/*
rm -r mod/*
rm run
make -f src/makefile.mak
rm *.o
rm -rf mod
echo "Make script complete..."
