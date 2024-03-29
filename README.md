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
./mkmc -k20 -thr_rat0.5 @<input_files> <output_file> <temp_dir>
```
It will generate a matrix of 20-mers occurring in at least a half of the input files.

#### General usage:
```
mkmc [options] @<input_files> <output_file> <temp_dir>
```

The main parameters are as follows:
 - `<input_files>` - file with a list of samples names with input files names in specified (-f<a/q/m> switch) format (gzipped or not)

File example:
```
killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz
killifishretina2 kfB.fastq.gz
```
 - `<output_file>` - file where the matrix of k-mers counts or FASTA file will be dumped
 - `<temp_dir>` - directory where temporary files will be stored		
 - `-thr<X>` - filter out k-mers occuring less than `<X>` times... (default: 1)
 - `-thr_rat<Y>` ... in a ratio `<Y>` of the input files (per k-mer sequence filtering) (default: 0.0)
 - `-flt<X>` - keep k-mers present in <X> file (FASTA or a sequence of the k-mers) only

E.g. `-thr_rat0.5` and `-thr2` mean that k-mers appearing at least twice in at least a half of the input files will be dumped. We recommend to be careful while specifying `<Y>` parameter, small values (e.g. 0.0) cause obtaining an enormous output file.

As `[options]` you can also pass optional parameters:
 - `-k<len>` - k-mer length (default: 25)
 - `-f<a/q/m>` - input in FASTA format (`-fa`), FASTQ format (`-fq`), or multi FASTA (`-fm`); mixing files is not supported (default: FASTQ)
 - `-of<a,matrix>` - output in FASTA format (`-ofa`) or matrix (`-ofmatrix`) (default: matrix)
 - `-b` - turn off transformation of k-mers into canonical form
 - `-ci<X>` - exclude counting k-mers occurring less than `<X>` times (if k-mer occurs less than `<X>` times in a file, it gets counter 0, but for this file only) (default: 1)
 - `-cx<X>` - exclude counting k-mers occurring more than `<X>` times (if k-mer occurs more than `<X>` times in a file, it gets counter 0, but for this file only) (default: 4e9)
 - `-cs<X>` - maximal value of a counter (default: 65535)
 - `-t<X>` - number of threads (default: 16 or no. of logic CPU cores)
 - `-wrk<X>` - number of parallel k-mer counting tasks (default: 4)
 - `-m<X>` - max amount of RAM in GB (from 2 to 1024); practically works only if `-r` is not set (default: 16)
 - `-r` - RAM only mode for k-mer counting
 - `-v` - verbose mode, shows progress

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
sample1 1.fq
sample2 2.fq
sample3 3.fq
```
All the samples are stored in single, unpaired files. To have k-mers that were present in each input sample one may use:
```
mkdir -p tmp
./mkmc -k25 -fq -thr_rat1 @files.txt present-in-all.txt tmp
```
The output (`present-in-all.txt`) is then:
```
k-mer	sample1	sample2	sample3	
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
k-mer	sample1	sample2	sample3	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```

To have k-mers that were present in any of the samples one may use the following:
```
mkdir -p tmp
./mkmc -k25 -fq -thr_rat0.0 @files.txt present-in-any.txt tmp
```
The output (`present-in-any.txt`) is then:
```
k-mer	sample1	sample2	sample3	
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
 - for all the 201 files ca. 3,5h of computation, 66GB of RAM, and 800GB of HDD for temporary files required, the output file (for `-thr_rat0.5`) size was 23GB
 - for a subset of 36 files ca. 50min. of computation, 50GB of RAM required
