#!/bin/bash

mkdir -p tmp
../bin/mkmc -k25 -fq -thr1 @files.txt present-in-all.txt tmp
../bin/mkmc -k25 -fq -thr0.5 @files.txt present-in-at-least-half-files.txt tmp
../bin/mkmc -k25 -fq -thr0 @files.txt present-in-any.txt tmp
