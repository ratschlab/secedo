source=$1 # contains full path to the source BAM file to be split
filename=$(basename -- "${source}")  # extract file name from path
pref="${filename%.*}"  # remove extension
base_dir=$2

for i in {1..22} X Y; do
  chunk_name="${base_dir}/aligned_cells_split/${pref}_chr${i}.bam"
  samtools view -b "${source}" ${i} > "${chunk_name}"
  samtools index "${chunk_name}"
done