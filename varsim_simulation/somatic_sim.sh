varsim_somatic.py --reference GRCh38_new.fa --id som --som_num_snp 40000 \
        --som_num_ins 250 --som_num_del 250 \
        --som_num_mnp 200 \
        --som_num_complex 200 \
        --cosmic_vcf cosmic.vcf.gz \
        --normal_vcf germ_out/germ.truth.vcf \
        --nlanes 700 --total_coverage 35 \
        --simulator_executable ../art_bin_MountRainier/art_illumina \
        --out_dir som_out --log_dir som_log --work_dir som_work \
        --sv_insert_seq empty_file &> somatic.log