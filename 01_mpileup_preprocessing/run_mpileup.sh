while read line; do
	command="samtools mpileup -BR -d 1000000 -r $line --output-QNAME -b bam_list | python3 preprocessing.py -n 2224 -t 0.001 > mpileups/breast_B_chr_"$line"_mpileup"
	echo bsub -W 120:00 "$command"
	bsub -W 120:00 "$command"
done < chr_names



