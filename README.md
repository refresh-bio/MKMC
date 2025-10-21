# MKMC (multi-KMC)
> [!warning]  
**!!! Currently, this repository contains work in progress and should probably not be used in production !!!**


MKMC is a software utilizing KMC to count k-mers in each of the predefined input samples.
Then it combines multiple KMC databases into one single matrix, which optionally is saved into binary .kmcdb file and a text file.
The latter file contains k-mers as rows and samples as columns. The values are counts of k-mers in samples.
FASTA output files, containg k-mers sequences only, are also supported.

The matrix may be utilized to compute many of statistics (normalization, correlation, differential k-mers analysis, cross-validation, entropy, UMAP, PCA).

## Building
The easiest way to get the program is to download the most recent version from the [**release page**](https://github.com/refresh-bio/MKMC/releases).

To build own binary clone the repository with the command:
```
git clone --recurse-submodules https://github.com/refresh-bio/MKMC-dev.git
```
To build on Linux type `make -j` (make and G++ 11 or newer are required). To build on Windows use Visual Studio 2022 or newer.

## Usage
```
./mkmc [OPTIONS] -- input_samples output_files_prefix temp_dir
```

Positionals:
  - `input_samples TEXT:FILE REQD` - file with a list of samples and input files in specified (`-f` parameter) format (gzipped or not)
  - `output_files_prefix TEXT REQD` - template (prefix) of output files names
  - `temp_dir TEXT:DIR REQD` - directory for temporary files

Options:
 - `-h,--help` - Print this help message and exit
 - `-k UINT:UINT in [1 - 256] [25]` - k-mer length
 - `--tot_cnt` - generate samples counts sums file
 
[Option Group: k-mers filtering]
  Options:
 - `--thr UINT:POSITIVE [1]` - filter out k-mers occuring less than specified number of times...
 - `--thr_rat FLOAT:FLOAT in [0 - 1] [0]` ... in a specified ratio of the input files (see example)
 - `--flt TEXT:FILE` - keep k-mers present in a specified file (FASTA or a set of the k-mers, one per line) only; `-b` is used accordingly
 
[Option Group: correlation and normalization]
  Options:
 - `-n ENUM:value in {deseq,freq,q}` - normalize counts (DESeq2/frequency count/quantile normalization) before use
 - `--save_n` - save matrix with normalized counts to file
 - `--cor ENUM:value in {kendall,pearson,spearman}` ... Needs: `-n` `-p` - compute correlation coefficients, basing on a phenotype file (Kendall Tau/Pearson/Spearman correlation)
 - `-p TEXT:FILE` Needs: `--cor` - set a phenotype file (a sequence of integers, one in each line)
 
[Option Group: differential k-mers analysis]
  Options:
 - `--diff ENUM:value in {anova,dids,snr,ttest,wrs}` ... Needs: `-c` - perform differential k-mers analysis (ANOVA, DIDS, Signal to Noise ratio, T-Test, Wilcoxon-rank sum (Mann-Whitney U test)); all except T-Test need `-n`; counts for T-Test are always unnormalized, increased by 1, and logarithmized
 - `--pval_corr ENUM:value in {b,bh,by,hb}` Needs: `--diff` - correct p-values of differential k-mers analysis (Bonferroni/Benjamini-Hochberg/Benjamini-Yekutieli/Holm-Bonferroni); store statistically significant k-mers also in separated files; useful for ANOVA, T-Test, Wilcoxon-rank sum
 - `--max_corrected_pval FLOAT:FLOAT in [0 - 1] [0.05]` Needs: `--pval_corr` - statistical significance for `--pval_corr` parameter
 - `-c TEXT:FILE` Needs: `--diff` - set a phenotype file for differential k-mers analysis (a sequence of natural numbers or text labels, one in each line)
 - `--dids-mode ENUM:value in {quadratic,sqrt,tanh} [sqrt]` - DIDS mode (x*x, square root, 1 + tanh(3x - 3))

[Option Group: cross-validation]
  Options:
 - `--cv` Needs: `--cor` - perform cross-validation for correlation
 - `--leave UINT [2]` Needs: `--cv` - number of samples to leave in every test
 - `--cv-seed UINT [1234567890]` Needs: `--cv` - random seed
 
[Option Group: other statistical parameters]
  Options:
 - `--entropy` - generate k-mers counts entropy; counts are increased by 1
 - `--n_top UINT [10000]` - select a maximal number of top k-mers by statistics with no p-values (for correlations in terms of an absolute value) and store them in separate files; needs `--cor` or `--diff`

[Option Group: dimentionality reduction]
  Options:
 - `--umap` Needs: `-n` - reduce dimensionality of normalized matrix with UMAP
 - `--pca` Needs: `-n` - reduce dimensionality of normalized matrix with PCA
 - `--dimensions UINT [2]` Needs: `--umap` or `--pca` - number of output dimensions
 - `--umap-local_connectivity FLOAT [1]` Needs: `--umap` - local_connectivity parameter
 - `--umap-bandwidth FLOAT [1]` Needs: `--umap` - `bandwidth` parameter
 - `--umap-mix_ratio FLOAT [1]` Needs: `--umap` - `mix_ratio` parameter
 - `--umap-spread FLOAT [1]` Needs: `--umap` - `spread` parameter
 - `--umap-min_dist FLOAT [0.01]` Needs: `--umap` - `min_dist` parameter
 - `--umap-a FLOAT [0]` Needs: `--umap` - `a` parameter
 - `--umap-b FLOAT [0]` Needs: `--umap` - `b` parameter
 - `--umap-repulsion_strength FLOAT [1]` Needs: `--umap` - `repulsion_strength` parameter
 - `--umap-initialize ENUM:value in {none,random,spectral,spectral_only} [spectral]` Needs: `--umap` - `initialize` parameter
 - `--umap-num_neighbors INT [15]` Needs: `--umap` - `num_neighbors` parameter
 - `--umap-num_epochs INT [-1]` Needs: `--umap` - `num_epochs` parameter
 - `--umap-learning_rate FLOAT [1]` Needs: `--umap` - `learning_rate` parameter
 - `--umap-negative_sample_rate FLOAT [5]` Needs: `--umap` - `negative_sample_rate` parameter
 - `--umap-seed UINT [1234567890]` Needs: `--umap` - `seed` parameter
 - `--umap-parallel_optimization INT [0]` Needs: `--umap` - `parallel_optimization` parameter
 - `--pca-mode ENUM:value in {covariance,svd} [svd]` Needs: `--pca` - PCA mode

[Option Group: additional parameters]
  Options:
 - `-f ENUM:value in {fa,fq,mf} [fq]` - input format (FASTA, FASTQ or multi-FASTA); mixing files formats is not supported
 - `-o ENUM:value in {fa,matrix} [matrix]  ...` - save k-mers (FASTA or matrix with unnormalized counts) to file
 - `-b` - turn off transformation of k-mers into canonical form; applies both for input sequences and k-mers passed by `--flt`
 - `--ci UINT:POSITIVE [1]` - exclude k-mers occurring less than specified number of times (if k-mer occurs less than --ci times in a sample, it gets counter 0, but for this sample only)
 - `--cx UINT:POSITIVE [4000000000]` - exclude counting k-mers occurring more than specified number of times (if k-mer occurs more than --cx times in a sample, it gets counter 0, but for this sample only)
 - `--cs UINT:UINT in [2 - 4294967295] [65535]` - maximal value of a counter
 - `--wrk UINT [4]` - number of parallel k-mer counting tasks
 - `-t UINT [no. of logic CPU cores]` - number of threads
 - `-m UINT:INT in [2 - 1024] [16]` - max amount of RAM in GB; practically works only if `-r` is not set
 - `-r` - count k-mers in RAM only
 - `-v` - verbose mode, shows progress and minor warnings, may be given up to 2 times

[Option Group: debug parameters]
  Options:
 - `--keep` - keep temporary files and binary results file
 - `--reuse-db` - reuse samples and filtering databases (if possible)
 - `--learn-deseq` Needs: `--keep` - collect data for DESeq2 normalization (not necessary for `-n deseq`, but useful for further  `--reuse-db`)
 - `--on UINT:POSITIVE [512]` - suggested number of internal bins, modify carefully
 - `--generate_snr_for_unnormalized_data` - generate Signal to Noise ratio also for unnormalized counts
K-mers order in output files is not specified and may vary between runnings.
> [!warning]  
**K-mers order in output files is not specified and may vary between runnings.**

## Examples
To obtain the statistics use i.a. one or many of the parameters: `-n`, `--cor`, `--diff`, `--cv`, `--entropy`, `--umap`, `--pca`. Statistics will be computed basing on counts matrix, which may be generated as follows.

### Generating matrix examples
```
./mkmc -k 20 --thr_rat 0.5 -- input_files_list.txt output tmp
```
It will generate a matrix (if `--keep` given, stored in a binary file `output.kmcdb`) of 20-mers occurring in at least a half of the input files.

```
./mkmc -k 20 --thr 2 --thr_rat 0.5 -- input_files_list.txt output tmp
```
It will generate a matrix of 20-mers occurring at least twice in at least a half of the input files.

To save the matrix to a text file use `-o matrix`:
```
./mkmc -k 20 --thr_rat 0.5 -o matrix -- input_files_list.txt output tmp
```
It will generate a file `output_matrix`.

`input_files_list.txt` example:
```
killifishretina1 kfA_1.fastq.gz kfA_2.fastq.gz
killifishretina2 kfB.fastq.gz`
```

### Results examples
The matrix is generated accordingly to the following. Let's assume FASTQ files contents:
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
The output matrix (refer `Generating matrix examples` section to see, how to obtain a text matrix file) is then:
```
k-mer	sample1	sample2	sample3	
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
```
To have k-mers that were present in at least half of the samples one may use the following:
```
./mkmc -k 25 -f fq --thr_rat 0.5 -- files.txt present-in-at-least-half-files tmp
```
The output matrix is then:
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
The output matrix is then:
```
k-mer	sample1	sample2	sample3	
AAAACACACAAACAGATAAACAGAT	1	1	0
ACCCCTGGGAATTAACCCACCTACG	0	0	1
ACCCCTGGGTTTTAACCCACGTACG	1	1	1
ACGTACGTGGGTTAAAACCCAGGGG	1	1	1
ACGTAGGTGGGTTAATTCCCAGGGG	0	0	1
TAAAACACACAAACAGATAAACAGA	1	1	0
```
