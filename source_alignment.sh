#Create alignment index
gmap_build -D ./genome -d hg38 -C hg38.fa

#Create splicesites file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
zcat gencode.v26.annotation.gtf.gz | gtf_splicesites | iit_store -o gencode.v26.splicesites.iit

#For each BAM file
for sample in `cat samples.txt`; do

  echo $sample;

  #Map sample against the reference genome
  gsnap -t 20 -B 4 -N 1 -E 1 -w 200000 --pairmax-rna 200000 -D genome -d hg38 -s gencode.v26.splicesites.iit $sample'_R1.fastq' $sample'_R2.fastq' --input-buffer-size 100000 --output-buffer-size 100000 --format sam 2>log_mapped_$sample'.txt' | samtools view -bS - > mapped_$sample'.bam';

  #Select high quality 101bp mapped reads (Q>20) and sort BAM file by genomic coordinates
  samtools view -q20 mapped_$sample.bam | awk '{if(length($10)==101) print$0}' | samtools view -bT genome/hg38.fa - | samtools sort - q20_mapped_$sample.sorted;

done
