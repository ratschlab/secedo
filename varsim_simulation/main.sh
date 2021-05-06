
module load jdk # varsim doesn't work with openjdk
module load bowtie2

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
coverage=0.01  # coverage for each generated cell
cov="cov${coverage#*.}x"  # e.g. cov01x
n_cells=500 # number of healthy cells and tumor cells in each group
n_tumor=4 # how many tumor cell types to generate

code_dir="$HOME/somatic_variant_calling/code"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

source ${SCRIPT_DIR}/1_run_generate_cell_pattern.sh
source ${SCRIPT_DIR}/2_run_generate_reads.sh
source ${SCRIPT_DIR}/3_run_mapping.sh
source ${SCRIPT_DIR}/4_run_pileup.sh
source ${SCRIPT_DIR}/5_run_variant_calling.sh


action=$1

if (( action == 1)); then
  generate_cell_patterns
fi
if (( action <= 2)); then
  generate_reads
fi
if (( action <= 3)); then
  map_reads
fi
if (( action <= 4)); then
  create_pileup
fi
if (( action <= 5)); then
  variant_calling
fi
