# projects2020-scvariant_calling
Master thesis for hanap

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
* -n <int> - number of cells in the data set
* -t <float> - theta, expected error rate

#### run_mpileup.sh
An example script to run samtools mpileup + preprocessing on each chromosome separately.

## 02_compute_simMat.py
Computes mat_same and mat_diff, matrices giving on position (i,j) the probability of observing reads from cells i and j assuming cells i and j have the same or
different genotypes, resp.

Required arguments:
* -n <int> - number of cells in the data set
* -e <float> - epsilon, estimated frequency of mutated loci in the pre-processed data set
* -h <float> - estimated frequency of homozygous loci in the pre-processed data set
* -i <string> - identifier of the run; the resulting files will be called "mat_same_<id>.csv" and "mat_diff_<id>.csv"
* -f <string> - input file, as outputted by preprocessing.py

Optional arguments:
* -t <float> - theta, expected error rate; default: 0.001
* -c <string> - file with identifiers of cells (numbers between 0 and (num_of_cells - 1); the matrices will be computed only for these cells; if not given, all cells are used
* -g <int> - file with cell groups: for each cell, an identifier of a group to which it belongs; cells from one group will be treated as one cell; useful e.g. when creating data with artificially higher coverage; if not given, each cell has its own group
* -m <int> - maximal considered insert size (for paired-end sequencing); default: 1000
