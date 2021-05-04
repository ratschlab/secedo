name=$1
pref="${name%.*}"
base_dir=$2

for i in {1..22} X Y; do
  chunk_name="${base_dir}/aligned_cells_split/${pref}_chr${i}.bam"
  samtools view -b "${base_dir}/aligned_cells/${name}" ${i} > "${chunk_name}"
  samtools index "${chunk_name}"
done