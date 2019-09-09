#Alternative Splicing Analysis using rMATS
python ./RNASeq-MATS.py -b1 q20_mapped_MES_83.sorted.bam,q20_mapped_MES_326.sorted.bam,q20_mapped_MES_1123.sorted.bam -b2 q20_mapped_PRO_19.sorted.bam,q20_mapped_PRO_157.sorted.bam,q20_mapped_PRO_528.sorted.bam -gtf gencode.v26.annotation.gtf -o results_rmats -t paired -analysis U -c 0.0001 -len 101
