# runs art_illumina to generate simulated reads for #n_cells healthy and #n_cells tumor cells with coverage #coverage

source global_vars.sh

# Start jobs for generating both healthy and tumor reads and wait for the jobs to complete
function generate_reads() {
  echo "Generating reads..."

  step=100

  art_illumina="/cluster/work/grlab/projects/projects2019-supervario/art_bin_MountRainier/art_illumina"
  gen_reads="python3 ~/somatic_variant_calling/code/varsim_simulation/generate_reads.py"
  scratch_dir="/scratch/svc"

  out_dir="${base_dir}/${cov}/healthy"
  mkdir -p "${out_dir}/logs/"
  out_prefix=${out_dir}/healthy_
  fasta="healthy.fa"

  for batch in $(seq 0 ${step} $((n_cells-1))); do
    cmd="echo Copying data...; mkdir -p ${scratch_dir}; cp ${base_dir}/genomes/${fasta}  ${scratch_dir}"
    cmd="$cmd;${gen_reads} --fasta ${scratch_dir}/${fasta} --art ${art_illumina} -p 20 --id_prefix healthy \
        --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage} \
        2>&1 | tee ${out_dir}/logs/sim-healthy-${batch}.log"
    echo ${cmd}
    bsub  -K -J "sim-he-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                -oo "${out_dir}/logs/sim-healthy-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}" &
  done

  out_dir="${base_dir}/${cov}/tumor"
  mkdir -p "${out_dir}/logs/"

  for tumor_type in $(seq 2 ${n_tumor}); do  # TODO: change back to 1
    out_prefix=${out_dir}/tumor_${tumor_type}_
    fasta="tumor-${tumor_type}.fa"

    for batch in $(seq 0 ${step} $((n_cells-1))); do
      cmd="echo Copying data...; mkdir -p ${scratch_dir}; cp ${base_dir}/genomes/${fasta} ${scratch_dir}"
      cmd="$cmd;${gen_reads} --fasta ${scratch_dir}/${fasta} --art ${art_illumina} -p 20 --id_prefix tumor_${tumor_type}_ \
          --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage} --seed_offset $((tumor_type*10000))  \
          2>&1 | tee ${out_dir}/logs/sim-tumor-${tumor_type}-${batch}.log"
      echo ${cmd}
      bsub  -K -J "sim-tu-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                  -oo "${out_dir}/logs/sim-tumor-${tumor_type}-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}" &
    done
  done

  wait

}