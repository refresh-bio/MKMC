# MKMC (multi-KMC)
> [!warning]  
**!!! Currently, this repository contains work in progress and should probably not be used in production !!!**


MKMC is a software utilizing KMC to count k-mers in each of the predefined input samples.
Then it combines multiple KMC databases into one single text file (currently, but probably will change in the future, i.e., more files, binary format).
This file is a matrix with k-mers as rows and samples as columns. Values are counts of k-mers in samples.
FASTA output files, containg k-mers sequences only, are also supported.

The easiest way to get the program is to download the most recent version from the [**release page**](https://github.com/refresh-bio/MKMC/releases).

#### General usage:
```
./mkmc [OPTIONS] input_samples_file output_files temp_dir
```

Positionals:
  - `input_samples_file TEXT:FILE REQUIRED` - file with a list of samples names with input files names in specified (`-f` parameter) format (gzipped or not)
  - `output_files TEXT REQUIRED` - file where the matrix of k-mers counts or FASTA file will be dumped
  - `temp_dir TEXT:DIR REQUIRED` - a directory where temporary files will be stored

Options:
 - `-h,--help` - Print this help message and exit
 - `-k UINT:UINT in [1 - 256] [25]` - k-mer length
 - `--thr UINT:POSITIVE [1]` -  filter out k-mers occuring less than specified number of times...
 - `--thr_rat FLOAT:FLOAT in [0 - 1] [0]` ... in a specified ratio of the input files (see example)
 - `--flt TEXT:FILE` - keep k-mers present in a specified file (FASTA or a set of the k-mers, one in each line) only
 - `-n ENUM:value in {freq,q}` - generate normalized counts (frequency count/quantile normalization)
 - `--cor ENUM:value in {kendall,pearson,spearman}` Needs: `-n` `-p` - compute correlation cofficients with specified methods, basing on a phenotype file (Kendall Tau/Pearson/Spearman correlation)	
 - `-p TEXT:FILE` Needs: `--cor` - set a phenotype file (a set of the integers, one in each line)

[Option Group: optional parameters]
  Options:	
 - `-f ENUM:value in {fa,fq,mf} [fq]` - input format (FASTA, FASTQ or multi-FASTA); mixing files is not supported
 - `-o ENUM:value in {fa,matrix} [matrix]  ...` - output format (FASTA or matrix)
 - `--on UINT:POSITIVE [512]` - number of output files, reduce carefully
 - `-b` - turn off transformation of k-mers into canonical form
 - `--ci UINT:POSITIVE [1]` - exclude k-mers occurring less than specified number of times (if k-mer occurs less than --ci times in a sample, it gets counter 0, but for this sample only)
 - `--cx UINT:POSITIVE [4000000000]` - exclude counting k-mers occurring more than specified number of times (if k-mer occurs more than --cx times in a sample, it gets counter 0, but for this sample only)
 - `--cs UINT:UINT in [2 - 4294967295] [65535]` - maximal value of a counter
 - `--wrk UINT [4]` - number of parallel k-mer counting tasks
 - `-t UINT [no. of logic CPU cores]` - number of threads
 - `-m UINT:INT in [2 - 1024] [16]` - max amount of RAM in GB; practically works only if `-r` is not set
 - `-r` - RAM only mode for k-mer counting
 - `-v` - verbose mode, shows progress

[Option Group: debug parameters]
  Options:
 - `--keep` - keep temporary files

Example: to run MKMC, type:
```
./mkmc -k 20 --thr_rat 0.5 input_files_list.txt output tmp
```
It will generate a matrix of 20-mers occurring in at least a half of the input files.
```
./mkmc -k 20 --thr 2 --thr_rat 0.5 input_files_list.txt output tmp
```
It will generate a matrix of 20-mers occurring at least twice in at least a half of the input files.

`input_files_list.txt` example:
```
killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz
killifishretina2 kfB.fastq.gz`
```

#### Results example
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
./mkmc -k25 -f fq -thr_rat 1 files.txt present-in-all.txt tmp
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
./mkmc -k 25 -f fq -thr_rat 0.5 files.txt present-in-at-least-half-files.txt tmp
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
./mkmc -k 25 -f fq -thr_rat 0.0 files.txt present-in-any.txt tmp
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
