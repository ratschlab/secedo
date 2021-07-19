DIR="/cluster/home/ddanciu/silver/datasets/gastric"
NEW_BAM=$DIR"/processed_files/mkn45_5k_750k_rpc_possorted_bam_filtered.bam"
NEW_BAM2=$DIR"/processed_files/mkn45_5k_750k_rpc_possorted_bam_filtered_withCBtag.bam"
HEADER=$DIR"/processed_files/tmp_header"
TMP=$DIR"/processed_files/tmp"
LOG="log_2_filteringCB"
echo > $LOG


#### filtering out reads not containing the CB tag
echo "Filtering reads without CB tag" >> $LOG
# save the header
samtools view -H $NEW_BAM > $HEADER
# filter the body
echo "samtools view $NEW_BAM | grep CB:Z: > $TMP" >> $LOG
samtools view -@ 4 $NEW_BAM | grep "CB:Z:"  > $TMP
echo "cat $HEADER $TMP | samtools view -b > $NEW_BAM2" >> $LOG
cat $HEADER $TMP | samtools view -b > $NEW_BAM2
echo "Number of reads with CB tag:" >> $LOG
samtools view -@ 4 -c $NEW_BAM2 >> $LOG
echo >> $LOG
# remove the temporary files
rm $TMP $HEADER
# index the file
echo "samtools index -@ 4 $NEW_BAM2" >> $LOG
samtools index -@ 4 $NEW_BAM2
