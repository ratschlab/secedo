DIR="/Volumes/Elements"
NEW_BAM="${DIR}/possorted_bam_filtered.bam"
NEW_BAM2="${DIR}/possorted_bam_filtered_withCBtag.bam"
HEADER="${DIR}/tmp_header"
TMP="${DIR}/tmp"
LOG="log_2_filteringCB"
echo > $LOG


#### filtering out reads not containing the CB tag
echo "Filtering reads without CB tag" >> $LOG

echo "samtools view -@ 4 -d CB -b $NEW_BAM >> $LOG"
samtools view -@ 4 -d CB -b $NEW_BAM > ${NEW_BAM2}

# index the file
echo "samtools index -@ 4 $NEW_BAM2" >> $LOG
samtools index -@ 4 $NEW_BAM2

echo "Number of reads with CB tag:" >> $LOG
samtools idxstats ${NEW_BAM2} | awk -F '\t' '{s+=$3+$4}END{print s}'
# samtools view -@ 4 -c $NEW_BAM2 >> $LOG
echo >> $LOG
