# MKMC (multi-KMC)

**!!! Currently this repository contains work in progress and should probably not be used in production!!!**

MKMC is software utylizing KMC to count k-mers in each of predefined input samples.
Then it combines multiple KMC databases into one single text file (currently, but probably will change in the future, i.e. more files, binary format).
This file is a matrix with k-mers as rows and samples as columns. Values are counts of k-mers in samples.

The easiest way to get the program is to download the most recent version from the [**release page**](https://github.com/refresh-bio/MKMC/releases).

To run MKMC on Linux (`mkmc` and `kmc_tools` have to be in the same directory) type:
```
./mkmc -k20 -thr0.5 [other_parameters] @<input_files> <output_file> <temp_dir>
```
It will generate a matrix of 20-mers from the input files.

#### The parameters are as follows

 - `<input_files>` is a text file with a list of FASTQ/FASTA files (one in each line)
 - `<output_file>` is a file where the matrix of k-mers counts will be dumped
 - `<temp_dir>` is a directory where temporary files will be stored
 - `-thr<X>` where <X> is a number of <0,1>, is the fraction of the input files a k-mer should be present in to be dumped into <output_file>; e.g. -thr0.5 mean, that k-mers appearing in at least a half of the input files will be dumped. We recommend to be careful while specifying this parameter, small values (e.g. 0) cause obtaining an enormous output file

As `[other_parameters]` you can also pass optional parameters:

 - `-fq` - if you would like to process FASTQ files (by default mkmc exptects FASTAs); mixing files types is not supported
 - `-wrk<X>` - in a case of a memory requirements reduction need (at the expense of a computation time) pass <X> smaller than 8, e.g. -wrk4

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
./mkmc -k25 -fq -thr1 @files.txt present-in-all.txt tmp
```
The output (`present-in-all.txt`) is then:
```
k-mer	1.fq	2.fq	3.fq	
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
```
To have k-mers that were present in at least half of samples one may use:
```
mkdir -p tmp
./mkmc -k25 -fq -thr0.5 @files.txt present-in-at-least-half-files.txt tmp
```
The output (`present-in-at-least-half-files.txt`) is then:
```
k-mer	1.fq	2.fq	3.fq	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```

To have k-mers that were present in any of samples one may use:
```
mkdir -p tmp
./mkmc -k25 -fq -thr0 @files.txt present-in-any.txt tmp
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

We have obtained the following benchmarks (on a computer equipped with AMD Ryzen Threadripper 3990X 64-Core processors):
 - for all the 201 files ca. 3,5h of computation, 66GB of RAM and 800GB of HDD for temporary files required, the output file (for -wrk0.5) size was 23GB
 - for a subset of 36 files ca. 50min. of computation, 50GB of RAM required
