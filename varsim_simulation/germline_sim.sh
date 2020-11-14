varsim.py --id germ --vc_in_vcf common_all_20180418.vcf.gz \
--reference GRCh38_new.fa \
--read_length 100 --vc_num_snp 3000000 --vc_num_ins 100000 \
--vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 \
--sv_num_ins 0 --sv_num_del 0 --sv_num_dup 0 --sv_num_inv 0 --sv_insert_seq empty_file \
--sv_dgv empty_file \
--mean_fragment_size 350 --sd_fragment_size 50 \
--nlanes 700 --total_coverage 35 --simulator_executable ../art_bin_MountRainier/art_illumina \
--out_dir germ_out --log_dir germ_log --work_dir germ_work