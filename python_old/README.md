# This is the original Python implementation of SVC and it doesn't reflect the C++ version anymore. We decided not to delete the code, as some Python code snippets may 
still be useful for quick prototyping of various ideas.

## What is where
* 00_filtering, ..., 06_varCalling - the individual steps of the algorithm (see below)
* varsim_simulation - the files used for simulating data using VarSim
* other - other useful scripts (see below)

## 00_filtering
Scripts used for filtering and formatting the 10x or simulated data for use with the rest of the pipeline.

Result: Each cell (defined by the "CB" tag for the 10x data) has its own sorted bam file. Only reads mapped in proper pair, non-duplicate and only primary alignments.

## 01_mpileup_preprocessing
The individual bam files are combined using samtools mpileup. Each output line is immediately pipelined to pre-processing (preprocessing.py) and saved to output
file only if it passes the pre-processing heuristics.

#### preprocessing.py
The code performing pre-processing.

Required arguments:
* -n \<int\> - number of cells in the data set
* -t \<float\> - theta, expected error rate

#### run_mpileup.sh
An example script to run samtools mpileup + preprocessing on each chromosome separately.

## 02_compute_simMat.py
Computes mat_same and mat_diff, matrices giving on position (i,j) the probability of observing reads from cells i and j assuming cells i and j have the same or
different genotypes, resp.

Required arguments:
* -n \<int\> - number of cells in the data set
* -e \<float\> - epsilon, estimated frequency of mutated loci in the pre-processed data set
* -h \<float\> - estimated frequency of homozygous loci in the pre-processed data set
* -i \<string\> - identifier of the run; the resulting files will be called "mat_same_\<id\>.csv" and "mat_diff_\<id\>.csv"
* -f \<string\> - input file, as outputted by preprocessing.py

Optional arguments:
* -t \<float\> - theta, expected error rate; default: 0.001
* -c \<string\> - file with identifiers of cells (numbers between 0 and (num_of_cells - 1); the matrices will be computed only for these cells; if not given, all cells are used
* -g \<int\> - file with cell groups: for each cell, an identifier of a group to which it belongs; cells from one group will be treated as one cell; useful e.g. when creating data with artificially higher coverage; if not given, each cell has its own group
* -m \<int\> - maximal considered insert size (for paired-end sequencing); default: 1000
  
## 03_compute_specClus.py
Combines the partial similarity matrices computed in the previous step, computes AIC and BIC to decide whether the data should be clustered, and clusters the data 
into two clusters using spectral clustering. The computed AIC and BIC and the labels reported by spectral clustering are saved in files "aic_bic_final" and "labels_specClus.csv", resp.

Required arguments:
* -n \<int\> - number of cells
* -i \<string\> - file containing comma separated identifiers of partial files; e.g. if the content of the file is "1,2", data will be read from matrices "mat_same_1.csv", "mat_diff_1.csv", "mat_same_2.csv" and "mat_diff_2.csv".
  
## 04_EM_stepA.py
One iteration of the EM algorithm - recomputes cluster priors, centers and likelihoods of belonging to each cluster for each cell. Typically, this would be run in parallel e.g. for each chromosome separately. 

Outputs two comma-separated text files, "loglikelihood_0_\<ind\>" and "loglikelihood_1_\<ind\>", where \<ind\> is the specified identifier, containing log-likelihood of belonging to cluster 0 and cluster 1, resp., for each cell.
  
  Required arguments:
  * -l \<string\> - comma-separated file containing soft labels of cells; labels are numbers between 0 and 1, giving the probability of belonging to cluster 1
  * -m \<string\> - the mpileup file, result of pre-processing
  * -i \<string\> - identifier to be used in names of output files
  * -t \<float\> - theta, expected error rate; default 0.001
  
## 05_EM_stepB.py
Combines together partial results from EM_stepA: Sums together log-likelihoods of being in cluster 0 and in cluster 1 over all parallel runs (e.g. over all chromosomes) and computes new cell labels.

Outputs a comma-separated file "new_labels" containing the new labels (=probabilities of being in cluster 1).

Required arguments:
* -l \<string\> - comma-separated file containing soft labels of cells; labels are numbers between 0 and 1, giving the probability of belonging to cluster 1
  
## 06_varCalling.py
Performs the variant calling.

Outputs a tab-separated file with information on the called variants.

Required arguments:
* -l \<string\> - comma-separated file containing soft labels of cells; labels are numbers between 0 and 1, giving the probability of belonging to cluster 1
* -f \<string\> - the mpileup file, result of pre-processing
* -i \<string\> - identifier, typically chromosome name; the output file will be called "varCalls_ch\<id\>"
  
Optional arguments:
* -t \<float\> - theta, error rate; default 0.01
* -h \<float\> - prior on heterozygous genotypes; default 0.001
  
## Other scripts
### create_new_pileup.py
In case of splitting to more than two clusters, this creates mpileup file containing data from only the specified cells. Pre-processing is run again for each line and some lines may be discarded.

Required arguments:
* -n \<string\> - name of the new mpileup file
* -o \<string\> - name of the original mpileup file
* -l \<string\> - name of file containing cell labels
* -L \<string\> - cells with this label will be included in the new file
* -c \<string\> - file containing identifiers of cells to be included in the new file
(either -l and -L, or -c are required)
  
 Optional arguments:
 * -t \<float\> - theta, error rate; default 0.001
  
 ### split_mpileup_file.py
 For some chromosomes, computation of the similarity matrix took too long. This script splits the mpileup file into smaller files. Each resulting file has at least 20000 lines and the split is made at such point, such that no read is split into two different files.
 
 ### varCalling_performance.py
 Compares the variant calls with true mutations - classifies each call as TP or FP, plus adds all FNs.
 
 Required arguments:
 * -i \<string\> - identifier, typically chromosome name; reads data from "varCalls_chr\<id\>", results are stored to "true_chr\<id\>"
