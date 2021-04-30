# germline cells
for i in $(seq 0 699); do
	OLD_BAM="cell_bams/cell"$i"_filtered_sorted.bam"
	HEADER="tmp_header"
	READ_NAMES="tmp_readNames"
	REST="tmp_rest"
	BODY="tmp_body"
	samtools view -H $OLD_BAM > $HEADER
	samtools view $OLD_BAM | cut -f 1 | sed 's/,//g' > $READ_NAMES
	samtools view $OLD_BAM | cut -f 1 --complement > $REST
	paste $READ_NAMES $REST > $BODY
	cat $HEADER $BODY | samtools view -b > $OLD_BAM
done 

# somatic cells
for i in $(seq 0 699); do
	OLD_BAM="cell_bams/som_cell"$i"_filtered_sorted.bam"
	HEADER="tmp_header"
	READ_NAMES="tmp_readNames"
	REST="tmp_rest"
	BODY="tmp_body"
	samtools view -H $OLD_BAM > $HEADER
	samtools view $OLD_BAM | cut -f 1 | sed 's/,//g' > $READ_NAMES
	samtools view $OLD_BAM | cut -f 1 --complement > $REST
	paste $READ_NAMES $REST > $BODY
	cat $HEADER $BODY | samtools view -b > $OLD_BAM
done 

rm $HEADER $READ_NAMES $REST $BODY