# runs art_illumina to generate simulated reads for #n_cells healthy and #n_cells tumor cells with coverage #coverage

module load jdk # varsim doesn't work with openjdk

cov="cov01x"
n_cells=500 # twice as many cells will be generated: healthy and tumor
step=100
coverage=0.05  # coverage for each generated cell

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
art_illumina="/cluster/work/grlab/projects/projects2019-supervario/art_bin_MountRainier/art_illumina"
gen_reads="python3 ~/somatic_variant_calling/code/varsim_simulation/generate_reads.py"
scratch_dir="/scratch/svc"

out_dir="${base_dir}/${cov}/healthy"
mkdir -p "${out_dir}/logs/"
out_prefix=${out_dir}/healthy_
fasta="${base_dir}/genomes/healthy.fa"

for batch in $(seq 0 ${step} ${n_cells}); do
  cmd="echo Copying data...; mkdir ${scratch_dir}; cp ${healthy_fa}  ${scratch_dir}"
  cmd="$cmd;${gen_reads} --fasta ${fasta} --art ${art_illumina} -p 20 \
  --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage}"
  echo ${cmd}
#  bsub  -J "sim-he-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
#              -oo "${out_dir}/logs/sim-healthy-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}"
done

out_dir="${base_dir}/${cov}/tumor"
mkdir -p "${out_dir}/logs/"
out_prefix=${out_dir}/tumor_
fasta="${base_dir}/genomes/tumor-1.fa"

for batch in $(seq 0 ${step} ${n_cells}); do
  cmd="echo Copying data...; mkdir ${scratch_dir}; cp ${tumor_fa} ${scratch_dir}"
  cmd="$cmd;${gen_reads} --fasta ${fasta} --art ${art_illumina} -p 20 \
  --start ${batch} --stop $((batch + step)) --out ${out_tumor} --coverage ${coverage}"
  echo ${cmd}
# bsub  -J "sim-tu-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
#              -oo "${out_dir}/logs/sim-tumor-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}"
done

