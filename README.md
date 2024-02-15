# MKMC (multi-KMC)
> [!warning]  
**!!! Currently, this repository contains work in progress and should probably not be used in production !!!**


MKMC is a software utilizing KMC to count k-mers in each of the predefined input samples.
Then it combines multiple KMC databases into one single text file (currently, but probably will change in the future, i.e., more files, binary format).
This file is a matrix with k-mers as rows and samples as columns. Values are counts of k-mers in samples.
FASTA output files, containg k-mers sequences only, are also supported.

The easiest way to get the program is to download the most recent version from the [**release page**](https://github.com/refresh-bio/MKMC/releases).

To run MKMC on Linux (`mkmc` and `kmc_tools` have to be in the same directory), type:
```
./mkmc -k20 -thr_rat0.5 [other_parameters] @<input_files> <output_file> <temp_dir>
```
It will generate a matrix of 20-mers from the input files.

#### The main parameters are as follows

 - `<input_files>` is a text file with a list of FASTQ/FASTA files (one in each line)
 - `<output_file>` is a file where the matrix of k-mers counts will be dumped
 - `<temp_dir>` is a directory where temporary files will be stored
 - `-thr_rat<X>` where `<X>` is a number of <0,1>, and
 - `-thr<Y>` where `<Y>` is natural number, is the fraction `<X>` of the input files a k-mer should be present at least `<Y>` times in to be dumped into <output_file>; e.g. -thr_rat0.5 and -thr2 mean, that k-mers appearing at least twice in at least a half of the input files will be dumped. We recommend to be careful while specifying `<X>` parameter, small values (e.g. 0) cause obtaining an enormous output file

As `[other_parameters]` you can also pass optional parameters:

 - `-ci<X>` - exclude counting k-mers occurring less than `<X>` times (if k-mer occurs less than <value> times in a file, it gets counter 0, but for this file only)
 - `-cx<X>` - exclude counting k-mers occurring more than `<X>` times (if k-mer occurs more than <value> times in a file, it gets counter 0, but for this file only)
 - `-fq` - if you would like to process FASTQ files (by default mkmc exptects FASTAs); mixing files types is not supported
 - `-ofa`/`-omatrix` - write output as a FASTA file or as a matrix (default)
 - `-r` - RAM only mode for k-mer counting
 - `-wrk<X>` -number of parallel k-mer counting tasks, currently should not significantly impact memory requirements (unless -r is set)
 - `-t<X>` - number of threads (default: no. of logic CPU cores)

#### Simple example
Lets assume following FASTQ files:
 - `1.fq`:
```
@Common k-mers
ACGTACGTGGGTTAAAACCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@k-mers in 1 and 2
ATCTGTTTATCTGTTTGTGTGTTTTA
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
 - `2.fq`:
```
@Common k-mers
ACGTACGTGGGTTAAAACCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@k-mers in 1 and 2
ATCTGTTTATCTGTTTGTGTGTTTTA
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
 - `3.fq`:
```
@Common k-mers
ACGTACGTGGGTTAAAACCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@k-mers only in 3
ACGTAGGTGGGTTAATTCCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```

And file `files.txt` containng:
```
1.fq
2.fq
3.fq
```
To have k-mers that were present in each input sample one may use:
```
mkdir -p tmp
./mkmc -k25 -fq -thr_rat1 @files.txt present-in-all.txt tmp
```
The output (`present-in-all.txt`) is then:
```
k-mer	1.fq	2.fq	3.fq	
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
```
To have k-mers that were present in at least half of the samples one may use the following:
```
mkdir -p tmp
./mkmc -k25 -fq -thr_rat0.5 @files.txt present-in-at-least-half-files.txt tmp
```
The output (`present-in-at-least-half-files.txt`) is then:
```
k-mer	1.fq	2.fq	3.fq	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```

To have k-mers that were present in any of the samples one may use the following:
```
mkdir -p tmp
./mkmc -k25 -fq -thr_rat0 @files.txt present-in-any.txt tmp
```
The output (`present-in-any.txt`) is then:
```
k-mer	1.fq	2.fq	3.fq	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGAATTAACCCACCTACG	0	0	1
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
ACGTAGGTGGGTTAATTCCCAGGGG	0	0	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```

#### Current performance
This is an initial version of code with limited optimizations and parallelism.

We have obtained the following benchmarks (on a computer equipped with AMD Ryzen Threadripper 3990X 64-Core processor) for 0.0.1 version:
 - for all the 201 files ca. 3,5h of computation, 66GB of RAM, and 800GB of HDD for temporary files required, the output file (for -wrk0.5) size was 23GB
 - for a subset of 36 files ca. 50min. of computation, 50GB of RAM required
