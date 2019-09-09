#For each BAM file
for sample in `cat samples.txt`; do

  echo $sample;

  #Run htseq to cluster mapped reads against the reference transcriptome
  htseq-count -f bam -m union -a 20 -r pos -s reverse q20_mapped_$sample'.sorted.bam' gencode.v26.annotation.gtf > counts_$sample'.txt' 2> log_counts_$sample'.txt';

done

#Merge count files
paste counts_MES_83.txt counts_MES_326.txt counts_MES_1123.txt counts_PRO_19.txt counts_PRO_157.txt counts_PRO_528.txt | cut -f1,2,4,6,8,10,12 | grep -v PAR_Y > counts_GSCs.txt

