module load jdk # varsim doesn't work with openjdk
module load bowtie2

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
coverage=0.01  # coverage for each generated cell
cov="cov${coverage#*.}x"  # e.g. cov01x
n_cells=500 # twice as many cells will be generated: healthy and tumor
code_dir="$HOME/somatic_variant_calling/code"
healthy_prefix="healthy"  # read id prefix for healthy cells
tumor_prefix="tumor2_" # read id prefix for tumor cells
