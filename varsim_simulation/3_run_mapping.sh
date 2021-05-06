# align synthetic reads (tumor+healthy) generated e.g. by varsim to the reference genome,
# then filter and sort the resulting BAM file
# indexing is necessary because we are splitting the resulting aligned cells later by chromosome (samtools can only
# split indexed files).

source global_vars.sh

# Starts jobs for mapping reads against the GRCh38 human genome and waits for the jobs to complete
function map_reads() {
  echo "Mapping reads using bowtie2, sorting and splitting by chromosome..."
  step=10

  mkdir -p "${base_dir}/${cov}/aligned_cells"
  mkdir -p "${base_dir}/${cov}/aligned_cells_split"

  for idx in $(seq 0 ${step} $((n_cells-1))); do
    cmd="echo hello"
    for i in $(seq "${idx}" $((idx+step-1))); do
      suf=$(printf "%03d" $i)

      file1="${base_dir}/${cov}/healthy/healthy_${i}.1.fq.gz"
      file2="${base_dir}/${cov}/healthy/healthy_${i}.2.fq.gz"
      bam_file="${base_dir}/${cov}/aligned_cells/healthy_${suf}.bam"
      cmd="${cmd}; bowtie2 -p 20 -x ${base_dir}/genomes/ref_index/GRCh38 -1 $file1 -2 $file2 \
           | samtools view -h -b -f 0x2 -F 0x500 - \
           | samtools sort -@ 10 -o ${bam_file}; samtools index ${bam_file}; \
           ${code_dir}/varsim_simulation/split.sh ${bam_file} ${base_dir}/${cov}"
    done
    # echo "${cmd}"
    bsub -K -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${base_dir}/logs/bowtie-healthy-${i}.lsf.log" "${cmd}" &
  done

  for tumor_type in $(seq 2 "${n_tumor}"); do  # TODO: change back to 1
    for idx in $(seq 0 ${step} $((n_cells-1))); do
        cmd="echo hello"
        for i in $(seq "${idx}" $((idx+step-1))); do
          suf=$(printf "%03d" ${i})

          file1="${base_dir}/${cov}/tumor/tumor_${tumor_type}_${i}.1.fq.gz"
          file2="${base_dir}/${cov}/tumor/tumor_${tumor_type}_${i}.2.fq.gz"
          bam_file="${base_dir}/${cov}/aligned_cells/tumor_${tumor_type}_${suf}.bam"
          cmd="${cmd}; bowtie2 -p 20 -x ${base_dir}/genomes/ref_index/GRCh38 -1 $file1 -2 $file2 \
              | samtools view -h -b -f 0x2 -F 0x500 - \
              | samtools sort -@ 10 -o ${bam_file}; samtools index ${bam_file}; \
              ${code_dir}/varsim_simulation/split.sh ${bam_file} ${base_dir}/${cov}"
        done
        # echo "${cmd}"
        bsub -K -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${base_dir}/logs/bowtie-tumor-${i}.lsf.log" "${cmd}" &
    done
  done

  wait

}