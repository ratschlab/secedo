for i in $(seq 0 699); do
	# somatic cell
	OLD_BAM="cell_bams/som_cell"$i"_filtered.bam"
	NEW_BAM="cell_bams/som_cell"$i"_filtered_sorted.bam"
	samtools sort -@ 10 -o $NEW_BAM $OLD_BAM 

	# germline cell
	OLD_BAM="cell_bams/cell"$i"_filtered.bam"
	NEW_BAM="cell_bams/cell"$i"_filtered_sorted.bam"
	samtools sort -@ 10 -o $NEW_BAM $OLD_BAM 

done