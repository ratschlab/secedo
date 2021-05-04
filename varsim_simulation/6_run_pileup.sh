# generate pileup files from all cells
work_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
cov="cov01x"
out_dir="${work_dir}/${cov}/pileups"
pileup="~/somatic_variant_calling/code/build/preprocess"

mkdir -p ${out_dir}
for chromosome in {1..22} X; do # Y was not added - maybe it confuses things
        scratch_dir="/scratch/pileup_${chromosome}"
        source_files=${work_dir}/${cov}/aligned_cells_split/*_chr${chromosome}.bam*
        num_files=`ls -l ${source_files} | wc -l`
        echo "Found ${num_files} files for chromosome ${chromosome}"
        copy_command="echo Copying data...; mkdir ${scratch_dir}; cp ${source_files} ${scratch_dir}"
        command="echo Running pileup binary...; $pileup -i ${scratch_dir}/ -o ${out_dir}/chromosome --num_threads 20 \
                --log_level=trace --min_base_quality 13 --max_coverage 1000 --seq_error_rate 0.001 \
                --chromosomes ${chromosome} | tee ${out_dir}/logs/pileup-${chromosome}.log"
        echo "Copy command: ${copy_command}"
        echo "Pileup command: $command"
        # allocating 40G scratch space; for the 1400 simulated Varsim cells, chromosomes 1/2 (the longest) need ~22G
        bsub  -J "pile-${chromosome}" -W 01:00 -n 20 -R "rusage[mem=4000,scratch=2000]" -R "span[hosts=1]" \
              -oo "${out_dir}/logs/pileup-${chromosome}.lsf.log" "${copy_command}; ${command}; rm -rf ${scratch_dir}"
done

