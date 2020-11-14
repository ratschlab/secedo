for i in $(seq 0 699); do
	# somatic cell
	OLD_BAM="cell_bams/som_cell"$i".sam"
	NEW_BAM="cell_bams/som_cell"$i"_filtered.bam"
	samtools view -h -b -f 0x2 -F 0x500 $OLD_BAM > $NEW_BAM

	# germline cell
	OLD_BAM="cell_bams/cell"$i".sam"
	NEW_BAM="cell_bams/cell"$i"_filtered.bam"
	samtools view -h -b -f 0x2 -F 0x500 $OLD_BAM > $NEW_BAM
done