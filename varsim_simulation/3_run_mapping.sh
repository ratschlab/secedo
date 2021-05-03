# align synthetic reads (tumor+healthy) generated e.g. by varsim to the reference genome,
# then filter and sort the resulting BAM file
module load bowtie2

step=10
cov="cov01x"
n_cells=1000

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
mkdir -p ${base_dir}/${cov}/aligned_cells
for idx in `seq 0 ${step} ${n_cells}`; do
  cmd="echo hello"
  for i in `seq ${idx} $((idx+step-1))`; do
    suf=$(printf "%03d" $i)

    file1="${base_dir}/${cov}/healthy/lane${i}.read1.fq.gz"
    file2="${base_dir}/${cov}/healthy/lane${i}.read2.fq.gz"
    bam_file="${base_dir}/${cov}/aligned_cells/healthy_${suf}.bam"
    cmd="${cmd}; bowtie2 -p 20 -x ${base_dir}/ref_index/GRCh38 -1 $file1 -2 $file2 | samtools view -h -b -f 0x2 -F 0x500 - | samtools sort -@ 10 -o ${bam_file}; samtools index ${bam_file}"
  done
#  echo ${cmd}
  bsub  -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${base_dir}/logs/bowtie-healthy-${i}.lsf.log" "${cmd}"
done

for idx in `seq 0 ${step} ${n_cells}`; do
    cmd="echo hello"
    for i in `seq ${idx} $((idx+step-1))`; do
      suf=$(printf "%03d" ${i})

      file1="${base_dir}/${cov}/tumor/lane${i}.read1.fq.gz"
      file2="${base_dir}/${cov}/tumor/lane${i}.read2.fq.gz"
      bam_file="${base_dir}/${cov}/aligned_cells/tumor_${suf}.bam"
      cmd="${cmd}; bowtie2 -p 20 -x ${base_dir}/ref_index/GRCh38 -1 $file1 -2 $file2 | samtools view -h -b -f 0x2 -F 0x500 - | samtools sort -@ 10 -o ${bam_file}; samtools index ${bam_file}"
    done
#    echo ${cmd}
    bsub  -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${base_dir}/logs/bowtie-tumor-${i}.lsf.log" "${cmd}"
done

