# runs art_illumina to generate simulated reads for #n_cells healthy and #n_cells tumor cells with coverage #coverage

module load jdk # varsim doesn't work with openjdk

cov="cov01x"
n_cells=500 # twice as many cells will be generated: healthy and tumor
step=100
coverage=0.01  # coverage for each generated cell

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
art_illumina="/cluster/work/grlab/projects/projects2019-supervario/art_bin_MountRainier/art_illumina"
gen_reads="python3 ~/somatic_variant_calling/code/varsim_simulation/generate_reads.py"
scratch_dir="/scratch/svc"

out_dir="${base_dir}/${cov}/healthy"
mkdir -p "${out_dir}/logs/"
out_prefix=${out_dir}/healthy_
fasta="healthy.fa"

for batch in $(seq 0 ${step} $((n_cells-1))); do
  cmd="echo Copying data...; mkdir -p ${scratch_dir}; cp ${base_dir}/genomes/${fasta}  ${scratch_dir}"
  cmd="$cmd;${gen_reads} --fasta ${scratch_dir}/${fasta} --art ${art_illumina} -p 20 \
  --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage} \
    2>&1 | tee ${out_dir}/logs/sim-healthy-${batch}.log"
  echo ${cmd}
#  bsub  -J "sim-he-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
#              -oo "${out_dir}/logs/sim-healthy-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}"
done

out_dir="${base_dir}/${cov}/tumor"
mkdir -p "${out_dir}/logs/"
out_prefix=${out_dir}/tumor_
fasta="tumor-1.fa"

for batch in $(seq 0 ${step} $((n_cells-1))); do
  cmd="echo Copying data...; mkdir -p ${scratch_dir}; cp ${base_dir}/genomes/${fasta} ${scratch_dir}"
  cmd="$cmd;${gen_reads} --fasta ${scratch_dir}/${fasta} --art ${art_illumina} -p 20 --seed_offset 10000 \
  --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage} 2>&1 | tee \
    ${out_dir}/logs/sim-tumor-${batch}.log"
  echo ${cmd}
  bsub  -J "sim-tu-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
              -oo "${out_dir}/logs/sim-tumor-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}"
done

