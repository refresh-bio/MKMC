#!/bin/bash

mkdir -p tmp
../bin/mkmc -k25 -fq -thr_rat1 @files.txt present-in-all.txt tmp
../bin/mkmc -k25 -fq -thr_rat0.5 @files.txt present-in-at-least-half-files.txt tmp
../bin/mkmc -k25 -fq -thr_rat0 @files.txt present-in-any.txt tmp
