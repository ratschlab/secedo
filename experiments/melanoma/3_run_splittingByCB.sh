DIR="/cluster/work/grlab/projects/projects2019-supervario/10x_data/"
NEW_BAM2=$DIR"processed_files/colo829_G1_1k_possorted_bam_filtered_withCBtag.bam"
NEW_BAM2_SORTED=$DIR"processed_files/colo829_G1_1k_possorted_bam_filtered_withCBtag_sorted.bam"
CELL_BAMS_DIR=$DIR"processed_files/cell_bams/"
TAGS_FILE=$DIR"processed_files/allowedTags"
PER_CELL_SUMMARY_FILE=$DIR"colo829_G1_1k_per_cell_summary_metrics.csv"
BAM_LIST=$CELL_BAMS_DIR"list_of_bams"
LOG="log_3_splittingByCB"
echo > $LOG

##### split based on CB tag to many smaller sam files
# first sort the file based on CB tag
echo "Sort by CB tag" >> $LOG
echo "samtools sort -t CB -@ 10 $NEW_BAM2 > $NEW_BAM2_SORTED" >> $LOG
samtools sort -t CB -@ 10 $NEW_BAM2 > $NEW_BAM2_SORTED
echo >> $LOG

# allowed tags
cat $PER_CELL_SUMMARY_FILE | cut -d ',' -f 1 | tail -n +2 > $TAGS_FILE

# then create a separate file for each allowed tag
echo "Splitting based on CB tag" >> $LOG
mkdir $CELL_BAMS_DIR
echo "python3 ~/sc_clustering/split_by_CBtag.py -f $NEW_BAM2_SORTED -o $CELL_BAMS_DIR -t $TAGS_FILE" >> $LOG
python3 ~/sc_clustering/split_by_CBtag.py -f $NEW_BAM2_SORTED -o $CELL_BAMS_DIR -t $TAGS_FILE

