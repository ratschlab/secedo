# runs Varsim to simulate 700 reads for healthy cells based on the GRCh38_new.fa reference genome

# don't forget to run: "conda activate py2" before running this script

module load jdk # varsim doesn't work with openjdk

cov="cov01x"
n_cells=1000
coverage=10  # 1000 lanes and coverage 10 coresponds to actual coverage per cell of 10/1000=0.01

base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
out_dir=${base_dir}/${cov}/healthy

touch ${out_dir}/empty_file

cmd="time ${base_dir}/varsim-0.8.4/varsim.py --id healthy --vc_in_vcf ${base_dir}/common_all_20180418.vcf.gz \
   --reference ${base_dir}/GRCh38_new.fa \
   --read_length 100 --vc_num_snp 3000000 --vc_num_ins 100000 \
   --vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 \
   --sv_num_ins 0 --sv_num_del 0 --sv_num_dup 0 --sv_num_inv 0 --sv_insert_seq ${out_dir}/empty_file \
   --sv_dgv empty_file \
   --mean_fragment_size 350 --sd_fragment_size 50 \
   --nlanes ${n_cells} --total_coverage ${coverage} \
   --simulator_executable /cluster/work/grlab/projects/projects2019-supervario/art_bin_MountRainier/art_illumina \
   --out_dir ${out_dir} --log_dir ${out_dir}/logs/ --work_dir ${out_dir}/tmp | tee 2>&1 ${out_dir}/healthy.log"

echo "Executing: $cmd"

bsub  -J "sim-healthy" -W 23:00 -n 20 -R "rusage[mem=8000]" -R "span[hosts=1]"  -oo "${out_dir}/healthy.lsf.log" "${cmd}"

