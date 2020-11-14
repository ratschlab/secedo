for i in $(seq 0 699); do
	file1="som_out/lane"$i".read1.fq.gz"
	file2="som_out/lane"$i".read2.fq.gz"
	samFile="cell_bams/som_cell"$i".sam"
	echo bowtie2 -p 20 -x ref_index/GRCh38 -1 $file1 -2 $file2 -S $samFile >> log_mapping_som
	bowtie2 -p 20 -x ref_index/GRCh38 -1 $file1 -2 $file2 -S $samFile
	echo $i finished >> log_mapping_som
done