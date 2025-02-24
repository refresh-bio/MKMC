# MKMC (multi-KMC)
> [!warning]  
**!!! Currently, this repository contains work in progress and should probably not be used in production !!!**


MKMC is a software utilizing KMC to count k-mers in each of the predefined input samples.
Then it combines multiple KMC databases into one single file binary .kmcdb file and, optionally, a text matrix.
The latter file is a matrix with k-mers as rows and samples as columns. The values are counts of k-mers in samples.
FASTA output files, containg k-mers sequences only, are also supported.

### Building
The easiest way to get the program is to download the most recent version from the [**release page**](https://github.com/refresh-bio/MKMC/releases).

To build own binary clone the repository with the command:
```
git clone --recurse-submodules https://github.com/refresh-bio/MKMC-dev.git
```
To build on Linux type `make -j` (make and G++ 11 or newer are required). To build on Windows use Visual Studio 2022 or newer.

### General usage
```
./mkmc [OPTIONS] -- input_samples_file output_files_template temp_dir
```

Positionals:
  - `input_samples_file TEXT:FILE REQUIRED` - file with a list of samples names with input files names in specified (`-f` parameter) format (gzipped or not)
  - `output_files_template TEXT REQUIRED` - template (prefix) of output files names
  - `temp_dir TEXT:DIR REQUIRED` - a directory where temporary files will be stored

Options:
 - `-h,--help` - Print this help message and exit
 - `-k UINT:UINT in [1 - 256] [25]` - k-mer length
 - `--tot_cnt` - generate samples counts sums file
 
[Option Group: k-mers filtering]
  Options:
 - `--thr UINT:POSITIVE [1]` -  filter out k-mers occuring less than specified number of times...
 - `--thr_rat FLOAT:FLOAT in [0 - 1] [0]` ... in a specified ratio of the input files (see example)
 - `--flt TEXT:FILE` - keep k-mers present in a specified file (FASTA or a set of the k-mers, one in each line) only; if `-b` is not set, the k-mers are converted to canonical form
 
[Option Group: correlation and normalization]
  Options:
 - `-n ENUM:value in {deseq,freq,q}` - generate normalized counts (DESeq2/frequency count/quantile normalization)
 - `--cor ENUM:value in {kendall,pearson,spearman}` ... Needs: `-n` `-p` - compute correlation cofficients, basing on a phenotype file (Kendall Tau/Pearson/Spearman correlation)
 - `-p TEXT:FILE` Needs: `--cor` - set a phenotype file (a sequence of integers, one in each line)
 
[Option Group: differential k-mers analysis]
  Options:
 - `--diff ENUM:value in {anova,dids,snr,ttest,wrs}` ... Needs: `-c` - perform differential k-mers analysis (ANOVA, DIDS, Signal to Noise ratio, T-Test, Wilcoxon-rank sum (Mann-Whitney U test)); all except T-Test need `-n`; counts for T-Test are increased by 1 and logarithmized
 - `--pval_corr ENUM:value in {b,bh,by,hb}` Needs: `--diff` - correct p-values of differential k-mers analysis (Bonferroni, Benjamini-Hochberg, Benjamini-Yekutieli, Holm-Bonferroni); store statistically significant k-mers also in separated files
 - `--max_corrected_pval FLOAT:FLOAT in [0 - 1] [0.05]` Needs: `--pval_corr` - statistical significance for --pval_corr parameter
 - `-c TEXT:FILE` Needs: `--diff` - set a phenotype file for differential k-mers analysis (a sequence of natural numbers or text labels, one in each line)
 
[Option Group: other statistical parameters]
  Options:
 - `--entropy` - generate k-mers counts entropy; counts are increased by 1
 - `--n_top UINT [10000]` - select a number of top k-mers (for correlations using an absolute value)

[Option Group: dimentionality reduction]
  Options:
 - `--umap` Needs: `-n` - run dimentionality reduction on normalized matrix with UMAP
 - `--pca` Needs: `-n` - run dimentionality reduction on normalized matrix with PCA
 - `--dimensions UINT [2]` Needs: `--umap` or `--pca` - number of output dimensions
 - `--umap-local_connectivity FLOAT [1]` Needs: `--umap` - local_connectivity parameter
 - `--umap-bandwidth FLOAT [1]` Needs: `--umap` - `bandwidth` parameter
 - `--umap-mix_ratio FLOAT [1]` Needs: `--umap` - `mix_ratio` parameter
 - `--umap-spread FLOAT [1]` Needs: `--umap` - `spread` parameter
 - `--umap-min_dist FLOAT [0.01]` Needs: `--umap` - `min_dist` parameter
 - `--umap-a FLOAT [0]` Needs: `--umap` - `a` parameter
 - `--umap-b FLOAT [0]` Needs: `--umap` - `b` parameter
 - `--umap-repulsion_strength FLOAT [1]` Needs: `--umap` - `repulsion_strength` parameter
 - `--umap-initialize ENUM:value in {none,random,spectral,spectral_only}` ... Needs: `--umap` - `initialize` parameter
 - `--umap-num_epochs INT [-1]` Needs: `--umap` - `num_epochs` parameter
 - `--umap-learning_rate FLOAT [1]` Needs: `--umap` - `learning_rate` parameter
 - `--umap-negative_sample_rate FLOAT [5]` Needs: `--umap` - `negative_sample_rate` parameter
 - `--umap-seed UINT [1234567890]` Needs: `--umap` - `seed` parameter
 - `--umap-parallel_optimization INT [0]` Needs: `--umap` - `parallel_optimization` parameter

[Option Group: additional parameters]
  Options:
 - `-f ENUM:value in {fa,fq,mf} [fq]` - input format (FASTA, FASTQ or multi-FASTA); mixing files formats is not supported
 - `-o ENUM:value in {fa,matrix} [matrix]  ...` - output format (FASTA or matrix)
 - `-b` - turn off transformation of k-mers into canonical form; applies both for input sequences and k-mers passed by `--flt`
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
 - `--on UINT:POSITIVE [512]` - number of internal bins, modify carefully
K-mers order in output files is not specified and may vary between runnings.
> [!warning]  
**K-mers order in output files is not specified and may vary between runnings.**

### Example
To run MKMC, type:
```
./mkmc -k 20 --thr_rat 0.5 -- input_files_list.txt output tmp
```
It will generate a binary matrix file `output.kmcdb` of 20-mers occurring in at least a half of the input files.

To obtain also a text matrix dump, type:
```
./mkmc -k 20 --thr_rat 0.5 -o matrix -- input_files_list.txt output tmp
```
It will generate also a file `output_matrix`.

```
./mkmc -k 20 --thr 2 --thr_rat 0.5 -- input_files_list.txt output tmp
```
It will generate a matrix of 20-mers occurring at least twice in at least a half of the input files.

`input_files_list.txt` example:
```
killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz
killifishretina2 kfB.fastq.gz`
```

### Results example
Let's assume following FASTQ files:
 - `1_1.fq`:
```
@Common k-mers
ACGTACGTGGGTTAAAACCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
 - `1_2.fq`:
```
@k-mers in 1 and 2
ATCTGTTTATCTGTTTGTGTGTTTTA
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
 - `2_1.fq`:
```
@Common k-mers
ACGTACGTGGGTTAAAACCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
- `2_2.fq`:
```
@k-mers in 1 and 2
ATCTGTTTATCTGTTTGTGTGTTTTA
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
 - `3_1.fq`:
```
@k-mers only in 3
ACGTAGGTGGGTTAATTCCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```
- `3_2.fq`:
```
@Common k-mers
ACGTACGTGGGTTAAAACCCAGGGGT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
```

And file `files.txt` containing:
```
sample1 1_1.fq 1_2.fq
sample2 2_1.fq 2_2.fq
sample3 3_1.fq 3_2.fq
```
To have k-mers that were present in each input sample one may use:
```
./mkmc -k25 -f fq --thr_rat 1 -- files.txt present-in-all tmp
```
The output (binary `present-in-all.kmcdb`; see `Example` section to see, how to obtain also a text matrix dump) is then:
```
k-mer	sample1	sample2	sample3	
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
```
To have k-mers that were present in at least half of the samples one may use the following:
```
./mkmc -k 25 -f fq --thr_rat 0.5 -- files.txt present-in-at-least-half-files tmp
```
The output (`present-in-at-least-half-files.kmcdb`) is then:
```
k-mer	sample1	sample2	sample3	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```

To have k-mers that were present in any of the samples one may use the following:
```
./mkmc -k 25 -f fq --thr_rat 0 -- files.txt present-in-any tmp
```
The output (`present-in-any.kmcdb`) is then:
```
k-mer	sample1	sample2	sample3	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGAATTAACCCACCTACG	0	0	1
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
ACGTAGGTGGGTTAATTCCCAGGGG	0	0	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```

### Current performance
This is an initial version of code with limited optimizations and parallelism.

We have obtained the following benchmarks (on a computer equipped with AMD Ryzen Threadripper 3990X 64-Core processor) for 0.0.1 version:
 - for all the 201 files ca. 3,5h of computation, 66GB of RAM, and 800GB of HDD for temporary files required, the output file (for `-thr_rat0.5`) size was 23GB
 - for a subset of 36 files ca. 50min. of computation, 50GB of RAM required
