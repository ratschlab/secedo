# Simulates data using Varsim+dwgsim, aligns it, piles it up and runs variant calling on it

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
coverage=0.05  # read coverage for each cell
cov="cov${coverage#*.}x"  # e.g. cov01x for coverage 0.01x
n_cells=500 # number of healthy and tumor cells in each group
n_tumor=4 # how many tumor cell types to generate
code_dir="$HOME/somatic_variant_calling/code"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


# runs Varsim to generate a pattern for healthy and tumor cells based on the GRCh38_new.fa reference genome
# Inputs:
#  - the human reference genome in fasta format e.g. GRCh38_new.fa
#  - common variations on the reference genome, e.g. common_all_20180418.vcf.gz
#  - a VCF file to draw mutations for the cancer cells from, e.g. cosmic.vcf.gz
# Outputs:
#  - one fasta file for the healthy cell and one fasta for the tumor cell

# don't forget to run: "conda activate py2" before running this script

# note the use of --disable_sim; this means that no reads are generated using Varsim (because it's slow and error
# prone). Instead we generate the reads ourself directly using art_illumina in the next step
function generate_cell_patterns() {
  module load jdk # varsim doesn't work with openjdk
  out_dir=${base_dir}/genomes
  mkdir -p  ${out_dir}

  touch "${out_dir}/empty_file"

  cmd="time ${base_dir}/varsim-0.8.4/varsim.py --id healthy --vc_in_vcf ${base_dir}/common_all_20180418.vcf.gz \
     --reference ${base_dir}/genomes/GRCh38_new.fa \
     --read_length 100 --vc_num_snp 3000000 --vc_num_ins 100000 \
     --vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 \
     --sv_num_ins 0 --sv_num_del 0 --sv_num_dup 0 --sv_num_inv 0 --sv_insert_seq ${out_dir}/empty_file \
     --sv_dgv empty_file \
     --disable_sim \
     --simulator_executable doesnt_matter_we_are_not_simulating \
     --out_dir ${out_dir} --log_dir ${out_dir}/logs/ --work_dir ${out_dir}/tmp | tee 2>&1 ${out_dir}/healthy.log"

  echo "Executing: $cmd"

  # bsub  -K -J "sim-healthy" -W 1:00 -n 10 -R "rusage[mem=8000]" -R "span[hosts=1]"  -oo "${out_dir}/healthy.lsf.log" "${cmd}" &


  # We have to wait for the healthy.truth.vcf to be generated by the varsim simulation for healthy cells
  # (because the tumor cells are based on healthy cells + cosmic mutations)

  printf "\n\n\n"
  healthy_vcf="${base_dir}/genomes/healthy.truth.vcf"
  if [ -f ${healthy_vcf} ]; then
    echo "Using existing ${healthy_vcf}"
  else
    echo -n "Waiting for ${healthy_vcf} to be generated..."
    while [ ! -f ${healthy_vcf} ]; do sleep 10; echo -n .;  done;
    echo "done"
  fi
  # TODO: replace --reference with ${base_dir}/genomes/GRCh38_new.fa to generate tumor cells based on the healthy
  # cells. Here we are simulating tumor cells that have 40K SNPs in common but only differ among themselves with 20K
  for i in $(seq 2 "${n_tumor}"); do # TODO: set back to 1
    command="time python2 ${base_dir}/varsim-0.8.4/varsim_somatic.py \
            --reference ${base_dir}/genomes/tumor-20K-1.fa \
            --id tumor-20K-${i} \
            --seed ${i} \
            --som_num_snp 20000 \
            --som_num_ins 250 \
            --som_num_del 250 \
            --som_num_mnp 200 \
            --som_num_complex 200 \
            --cosmic_vcf ${base_dir}/cosmic/cosmic.vcf.gz \
            --normal_vcf ${healthy_vcf} \
            --disable_sim \
            --simulator_executable ${out_dir}/empty_file \
            --out_dir ${out_dir} \
            --log_dir ${out_dir}/logs --work_dir ${out_dir}/tmp --sv_insert_seq ${out_dir}/empty_file"

    echo "Executing: ${command}"

    # takes about 15 minutes
    bsub  -K -J "sim-tumor-${i}" -W 1:00 -n 10 -R "rusage[mem=20000]" -R "span[hosts=1]"  -oo "${out_dir}/logs/tumor-${i}.lsf.log" "${command}" &
  done

  wait
}

# runs art_illumina to generate simulated reads for #n_cells healthy and #n_cells tumor cells with coverage #coverage
# Start jobs for generating both healthy and tumor reads and wait for the jobs to complete
function generate_reads() {
  echo "Generating reads..."

  module load bowtie2

  step=100

  art_illumina="/cluster/work/grlab/projects/projects2019-supervario/art_bin_MountRainier/art_illumina"
  gen_reads="python3 ${code_dir}/varsim_simulation/generate_reads.py"
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
#    echo ${cmd}
    bsub  -K -J "sim-he-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                -oo "${out_dir}/logs/sim-healthy-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}" &
  done

  out_dir="${base_dir}/${cov}/tumor"
  mkdir -p "${out_dir}/logs/"

  for tumor_type in $(seq 1 "${n_tumor}"); do
    out_prefix=${out_dir}/tumor_${tumor_type}_
    fasta="tumor-20K-${tumor_type}.fa"

    for batch in $(seq 0 ${step} $((n_cells-1))); do
      cmd="echo Copying data...; mkdir -p ${scratch_dir}; cp ${base_dir}/genomes/${fasta} ${scratch_dir}"
      cmd="$cmd;${gen_reads} --fasta ${scratch_dir}/${fasta} --art ${art_illumina} -p 20 --id_prefix tumor_${tumor_type}_ \
          --start ${batch} --stop $((batch + step)) --out ${out_prefix} --coverage ${coverage} --seed_offset $((tumor_type*10000))  \
          2>&1 | tee ${out_dir}/logs/sim-tumor-${tumor_type}-${batch}.log"
 #     echo ${cmd}
      bsub  -K -J "sim-tu-${tumor_type}-${batch}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                  -oo "${out_dir}/logs/sim-tumor-${tumor_type}-${batch}.lsf.log" "${cmd}; rm -rf ${scratch_dir}" &
    done
  done

  wait
}


# align synthetic reads (tumor+healthy) generated e.g. by varsim to the reference genome,
# then filter and sort the resulting BAM file
# indexing is necessary because we are splitting the resulting aligned cells later by chromosome (samtools can only
# split indexed files).

# Starts jobs for mapping reads against the GRCh38 human genome and waits for the jobs to complete
# takes < 10 minutes
function map_reads() {
  echo "Mapping reads using bowtie2, sorting and splitting by chromosome..."
  step=10

  mkdir -p "${base_dir}/${cov}/aligned_cells"
  mkdir -p "${base_dir}/${cov}/aligned_cells_split"
  logs_dir="${base_dir}/${cov}/aligned_cells_split/logs"
  mkdir -p ${logs_dir}

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
    bsub -K -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${logs_dir}/bowtie-healthy-${i}.lsf.log" "${cmd}" &
  done

  for tumor_type in $(seq 1 "${n_tumor}"); do  # TODO: change back to 1
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
        bsub -K -J "bt-${i}" -W 2:00 -n 20 -R "rusage[mem=800]" -R "span[hosts=1]"  -oo "${logs_dir}/bowtie-tumor-${i}.lsf.log" "${cmd}" &
    done
  done

  wait

}

# Starts jobs for creating pileup files from the aligned BAM files. One job per Chromosome.
# Waits for jobs to complete
function create_pileup() {
  echo "Generating pileups..."
  out_dir="${base_dir}/${cov}/pileups"
  log_dir="${out_dir}/logs"
  pileup="${code_dir}/build/preprocess"

  mkdir -p ${out_dir}
  mkdir -p ${log_dir}
  for chromosome in {1..22} X; do # Y was not added - maybe it confuses things
          scratch_dir="/scratch/pileup_${chromosome}"
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
  echo "Running variant calling..."
  module load openblas
  work_dir="${base_dir}/${cov}"
  input_dir="${work_dir}/pileups"
  svc="${code_dir}/build/svc"
  flagfile="${code_dir}/flags_sim"
  out_dir="${work_dir}/svc/"
  mkdir -p "${out_dir}"
  command="${svc} -i ${input_dir}/ -o ${out_dir} --num_threads 20 --log_level=trace --flagfile ${flagfile} \
           --clustering_type SPECTRAL6 --merge_count 1 --max_coverage 100 | tee ${out_dir}/svc.log"
  echo "$command"

  bsub -K -J "svc" -W 03:00 -n 20 -R "rusage[mem=20000]" -R "span[hosts=1]" -oo "${out_dir}/svc.lsf.log" "${command}"
}

# check the command-line arguments
if [ "$#" -ne 1 ]; then
            echo_err "Usage: main.sh <start_step>"
            echo "start_step=1 -> Generate genomes for healthy/tumor cells (~30 mins)"
            echo "start_step=2 -> Generate reads for healthy/tumor cells (~20 mins)"
            echo "start_step=3 -> Align reads against the human genome (~10 mins)"
            echo "start_step=4 -> Create pileup files (one per chromosome) (~10 mins)"
            echo "start_step=5 -> Run variant calling (~10 mins)"
            exit 1
fi

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
