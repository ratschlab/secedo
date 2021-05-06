# generate pileup files from all cells

source global_vars.sh

# Starts jobs for creating pileup files from the aligned BAM files. One job per Chromosome.
# Waits for jobs to complete
function create_pileup() {
  echo "Generating pileups..."
  out_dir="${base_dir}/${cov}/pileups"
  log_dir="${out_dir}/logs"
  pileup="~/somatic_variant_calling/code/build/preprocess"

  mkdir -p ${out_dir}
  for chromosome in {1..22} X; do # Y was not added - maybe it confuses things
          scratch_dir="/scratch/pileup_${chromosome}"
          source_files=${base_dir}/${cov}/aligned_cells_split/*_chr${chromosome}.bam*
          num_files=`ls -l ${source_files} | wc -l`
          echo "Found ${num_files} files for chromosome ${chromosome}"
          copy_command="echo Copying data...; mkdir ${scratch_dir}; cp ${source_files} ${scratch_dir}"
          command="echo Running pileup binary...; $pileup -i ${scratch_dir}/ -o ${out_dir}/chromosome --num_threads 20 \
                  --log_level=trace --min_base_quality 13 --max_coverage 1000 --seq_error_rate 0.001 \
                  --chromosomes ${chromosome} | tee ${log_dir}/pileup-${chromosome}.log"
          echo "Copy command: ${copy_command}"
          echo "Pileup command: $command"
          # allocating 40G scratch space; for the 1400 simulated Varsim cells, chromosomes 1/2 (the longest) need ~22G
          bsub  -K -J "pile-${chromosome}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
                -oo "${log_dir}/pileup-${chromosome}.lsf.log" "${copy_command}; ${command}; rm -rf ${scratch_dir}" &
  done

  wait
}