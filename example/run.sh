#!/bin/bash

mkdir -p tmp
../bin/mkmc -k25 -f fq -thr_rat 1 files.txt present-in-all.txt tmp
../bin/mkmc -k 25 -f fq -thr_rat 0.5 files.txt present-in-at-least-half-files.txt tmp
../bin/mkmc -k 25 -f fq -thr_rat 0 files.txt present-in-any.txt tmp
