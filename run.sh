#!/bin/bash
echo "Run script started..."
mkdir sim/out
rm -r sim/out/*
rm -rf mod
rm *.o
./run
echo "Run script complete..."
