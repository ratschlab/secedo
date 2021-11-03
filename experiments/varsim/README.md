This directory contains the scripts to generate synthetic data using varsim and cluster that data.
## Conda environment
Varsim requires python 2. The following conda packages are needed:
  - pyvcf
  - numpy
  - pysam
  - pybedtools
```
conda create --name py2 python=2.7 pyvcf numpy pysam pybedtools
conda activate py2
```

Another weirdness of varsim is that it doesn't work with OpenJDK. You need the Oracle JDK.

## Creating genomes
Tweak the `create_genomes.sh` script to your needs. This script uses the GRCH38 reference genome, common variations from
dbSNP and the Cosmic database to create the healthy and tumor cells. All 3 of the above are available for download.

## Generate reads, align reads, create pileups and variant call
The `variant_call.sh` script does all of the steps above. Edit it and tweak it to your needs.

## Other scripts
 - The `generate_reads.py` script is called by `variant_call.sh` to generate reads using dwgsim based on fasta templates
for cancer and healthy cells
 - `split.sh` - simple script that splits a BAM by chromosome for distributed processing
 - `cov_vs_snp.sh` - used to test the lowest SNP count/coverage at which SECEDO is able to cluster 1000 cells (Used 
   to create Figure 5 in the paper).
