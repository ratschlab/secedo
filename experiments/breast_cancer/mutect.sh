# Calls Mutect on the BAM files corresponding to the clusters found by SILVER
slice="B"

silver_dir="/cluster/work/grlab/projects/projects2019-secedo"
bam_dir="${silver_dir}/datasets/breastcancer/slice${slice}/processed_files/cluster_bams"
mutect="~/jre1.6.0_45/bin/java -Xmx8g -jar ~/mutect/muTect-1.1.4.jar --analysis_type MuTect"
genome_dir="${silver_dir}/genomes"
log_dir="${bam_dir}/logs"

# check the command-line arguments
if [ "$#" -ne 1 ]; then
            echo "Usage: mutect.sh <cluster_number>"
            exit 1
fi

cluster=$1

mkdir -p ${log_dir}

for chromosome in {1..22} X Y; do

cmd="${mutect} \
  --reference_sequence ${genome_dir}/GRCh37.p13.genome.fa \
  --dbsnp ${genome_dir}/dbsnp_132_b37.leftAligned_1Y.vcf \
  --cosmic ${genome_dir}/cosmic_v94_hg37_coding_and_noncoding.vcf \
  --input_file:normal ${bam_dir}/clone19_1Y.bam \
  --input_file:tumor ${bam_dir}/clone${cluster}_1Y.bam \
  --out ${bam_dir}/clone${cluster}_${chromosome}.txt \
  --coverage_file ${bam_dir}/clone${cluster}_${chromosome}.coverage \
  --intervals $chromosome"

bsub -J "mutect_${cluster}_${chromosome}" -W 04:00 -n 2 -R "rusage[mem=8000]" \
  -R  "span[hosts=1]" -oo "${log_dir}/mutect_${cluster}_${chromosome}.lsf.log" "${cmd}"

echo ${cmd}

done
