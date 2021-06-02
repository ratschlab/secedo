# Simulates data using Varsim+dwgsim, aligns it, piles it up and runs variant calling on it

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
coverage=0.05  # read coverage for each cell
cov="cov${coverage#*.}x8K"  # e.g. cov01x for coverage 0.01x

code_dir="$HOME/somatic_variant_calling/code"

genomes=("${base_dir}/genomes/healthy.fa" "${base_dir}/genomes/tumor-40K-1/tumor-40K-1.fa" \
"${base_dir}/genomes/tumor-20K-2/tumor-20K-2.fa" "${base_dir}/genomes/tumor-20K-3/tumor-20K-3.fa" \
"${base_dir}/genomes/tumor-10K-4/tumor-10K-4.fa" "${base_dir}/genomes/tumor-15K-7/tumor-15K-7.fa" \
"${base_dir}/genomes/tumor-5K-5/tumor-5K-5.fa" "${base_dir}/genomes/tumor-2.5K-9/tumor-2.5K-9.fa")
n_tumor=${#genomes[@]}

n_cells=(2000 500 500 500 500 2000 500 1000) # number of cells in each group

# runs art_illumina to generate simulated reads for each group of cells with coverage #coverage
function generate_reads() {
  echo "[$(date)] Generating reads..."

  step=100

  art_illumina="/cluster/work/grlab/projects/projects2019-supervario/art_bin_MountRainier/art_illumina"
  gen_reads="python3 ${code_dir}/experiments/varsim/generate_reads.py"

  out_dir="${base_dir}/${cov}/cells"
  mkdir -p "${out_dir}/logs/"

  for cell_idx in $(seq 0 $(( n_tumor-1 ))); do
    scratch_dir=$(mktemp -d -t fasta-XXXXXXXXXX --tmpdir=/scratch)
    out_prefix=${out_dir}/cell_${cell_idx}_
    fasta="${genomes[${cell_idx}]}"
    fasta_fname=$(basename -- "${fasta}")  # extract file name from path

    for batch in $(seq 0 ${step} $((n_cells[cell_idx]-1))); do
      cmd="echo [$(date)] Copying data...; mkdir -p ${scratch_dir}; cp ${fasta} ${scratch_dir}"
      cmd="$cmd;${gen_reads} --fasta ${scratch_dir}/${fasta_fname} --art ${art_illumina} -p 20 --id_prefix \
      cell_${cell_idx}_  --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage} \
      --seed_offset $((cell_idx*10000))  2>&1 | tee ${out_dir}/logs/sim-tumor-${cell_idx}-${batch}.log"
      echo ${cmd}
      bsub  -K -J "sim-tu-${cell_idx}-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                  -oo "${out_dir}/logs/sim-tumor-${cell_idx}-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}" &
    done
  done

  wait
}


# align synthetic reads (tumor+healthy) generated e.g. by varsim to the reference genome,
# then filter and sort the resulting BAM file
# indexing is necessary because we are splitting the resulting aligned cells later by chromosome (samtools can only
# split indexed files).

# Starts jobs for mapping reads against the GRCh38 human genome and waits for the jobs to complete
# takes ~ 10 minutes for 1000 cells
function map_reads() {
  module load bowtie2
  echo "[$(date)] Mapping reads using bowtie2, sorting and splitting by chromosome..."
  step=10

  mkdir -p "${base_dir}/${cov}/aligned_cells"
  mkdir -p "${base_dir}/${cov}/aligned_cells_split"
  logs_dir="${base_dir}/${cov}/aligned_cells_split/logs"
  mkdir -p ${logs_dir}

  for cell_idx in $(seq 0 $(( n_tumor-1 ))); do
    for idx in $(seq 0 ${step} $((n_cells[cell_idx]-1))); do
        cmd="echo hello"
        for i in $(seq "${idx}" $((idx+step-1))); do
          suf=$(printf "%03d" ${i})

          file1="${base_dir}/${cov}/cells/cell_${cell_idx}_${i}.1.fq.gz"
          file2="${base_dir}/${cov}/cells/cell_${cell_idx}_${i}.2.fq.gz"
          bam_file="${base_dir}/${cov}/aligned_cells/cell_${cell_idx}_${suf}.bam"
          cmd="${cmd}; bowtie2 -p 20 -x ${base_dir}/genomes/ref_index/GRCh38 -1 $file1 -2 $file2 \
              | samtools view -h -b -f 0x2 -F 0x500 - \
              | samtools sort -@ 10 -o ${bam_file}; samtools index ${bam_file}; \
              ${code_dir}/experiments/varsim/split.sh ${bam_file} ${base_dir}/${cov}"
        done
        # echo "${cmd}"
        bsub -K -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${logs_dir}/bowtie-tumor-${i}.lsf.log" "${cmd}" &
    done
  done

  wait

}

# Starts jobs for creating pileup files from the aligned BAM files. One job per Chromosome.
# Waits for jobs to complete
function create_pileup() {
  module load openblas

  echo "[$(date)] Generating pileups..."
  out_dir="${base_dir}/${cov}/pileups"
  log_dir="${out_dir}/logs"
  pileup="${code_dir}/build/preprocess"

  mkdir -p ${out_dir}
  mkdir -p ${log_dir}
  for chromosome in {1..22} X; do # Y was not added - maybe it confuses things
          scratch_dir=$(mktemp -d -t pileup-XXXXXXXXXX --tmpdir=/scratch)
          source_files=${base_dir}/${cov}/aligned_cells_split/*_chr${chromosome}.bam*
          num_files=`ls -l ${source_files} | wc -l`
          echo "Found ${num_files} files for chromosome ${chromosome}"
          copy_command="echo Copying data...; mkdir ${scratch_dir}; cp ${source_files} ${scratch_dir}"
          command="echo Running pileup binary...; ${pileup} -i ${scratch_dir}/ -o ${out_dir}/chromosome --num_threads 20 \
                  --log_level=trace --min_base_quality 30 --max_coverage 1000 \
                  --chromosomes ${chromosome} | tee ${log_dir}/pileup-${chromosome}.log"
          echo "Copy command: ${copy_command}"
          echo "Pileup command: $command"
          # allocating 40G scratch space; for the 1400 simulated Varsim cells, chromosomes 1/2 (the longest) need ~22G
          bsub  -K -J "pile-${chromosome}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                -oo "${log_dir}/pileup-${chromosome}.lsf.log" "${copy_command}; ${command}; rm -rf ${scratch_dir}" &
  done

  wait
}

# Runs the variant caller on pileup files (either binary generated by our own pileup binary, or textual
# generated by samtools mpileup)
# ~5 minutes for 1000 cells coverage 0.05x
function variant_calling() {
  echo "[$(date)] Running variant calling..."
  module load openblas
  work_dir="${base_dir}/${cov}"
  input_dir="${work_dir}/pileups"
  svc="${code_dir}/build/svc"
  flagfile="${code_dir}/flags_sim"
  for hprob in 0.5; do
    for seq_error_rate in 0.05; do
      out_dir="${work_dir}/svc_${hprob#*.}_${seq_error_rate#*.}/"
      mkdir -p "${out_dir}"
      command="${svc} -i ${input_dir}/ -o ${out_dir} --num_threads 20 --log_level=trace --flagfile ${flagfile} \
             --homozygous_prob=${hprob} --seq_error_rate=${seq_error_rate} \
             --clustering_type SPECTRAL6 --merge_count 1 --max_coverage 1000 | tee ${out_dir}/svc.log"
             #       --pos_file=${base_dir}/cosmic/cosmic.vcf \
      echo "$command"

      bsub -K -J "svc" -W 08:00 -n 20 -R "rusage[mem=80000]" -R "span[hosts=1]" -oo "${out_dir}/svc.lsf.log" "${command}" &
    done
  done

  wait
}

# check the command-line arguments
if [ "$#" -ne 1 ]; then
            echo "Usage: main.sh <start_step>"
            echo "start_step=1 -> Generate reads for healthy/tumor cells (~20 mins)"
            echo "start_step=2 -> Align reads against the human genome (~10 mins)"
            echo "start_step=3 -> Create pileup files (one per chromosome) (~10 mins)"
            echo "start_step=4 -> Run variant calling (~10 mins)"
            exit 1
fi

action=$1

if (( action <= 1)); then
  generate_reads
fi
if (( action <= 2)); then
  map_reads
fi
if (( action <= 3)); then
  create_pileup
fi
if (( action <= 4)); then
  variant_calling
fi
