# TODO: remove this once we switch to our own mpileup and move indexing to run_mapping
# Removes commas from the read ids of the reads generated with ./varsim_somatic.py (i.e. tumor reads), as these
# commas are mishandled by samtools mpileup. Also indexes the file 
# Commas are added by varsim to separate regions in the read that map to the reference from deletions/insertions

cov="cov05x"
for f in ${cov}/aligned_cells/healthy_*.bam; do
    HEADER="${f}_header"
    READ_NAMES="${f}_readNames"
    REST="${f}_rest"
    BODY="${f}_body"

    cmd="samtools view -H ${f} > $HEADER; \
      samtools view ${f} | cut -f 1 | sed 's/,//g' > $READ_NAMES; \
      samtools view ${f} | cut -f 1 --complement > $REST; \
      paste $READ_NAMES $REST > $BODY; \
      cat $HEADER $BODY | samtools view -b > ${f}; \
      rm ${HEADER} ${READ_NAMES} ${BODY} ${REST};"

    bsub -J "${f}" -W 0:30 -n 1 -R "rusage[mem=8000]" -R "span[hosts=1]"  -oo "logs/comma-${f}.lsf.log" "${cmd}"
done

