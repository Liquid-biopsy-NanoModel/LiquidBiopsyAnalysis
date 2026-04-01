#!/bin/bash
cd /home/ssingh6/Seqdata/1612/
mkdir bedfiles
module load samtools
samtools index 1612.hac_hmcg.aligned.sort.bam
modkit pileup 1612.hac_hmcg.aligned.sort.bam ./bedfiles/1612_pileup.bed --ref /home/ssingh6/nanoDx/references/hg38.fa  --preset traditional --only-tabs  --edge-filter 0,27 --threads 20

#liftOver ./bedfiles/1954_pileup.bed /home/ssingh6/nanoDx/static/hg38ToHg19.over.chain.gz ./bedfiles/1954_pileup_hg19.bed unmapped.bed
cut -f1-11 ./bedfiles/1612_pileup.bed > ./bedfiles/temp_a.bed

bedtools intersect -a ./bedfiles/temp_a.bed -b /home/ssingh6/nanoDx/static/450K_hg38.bed -wa -wb | awk -v OFS='\t' '$4=$15' | cut -f1-11 > ./bedfiles/1612_CpG.450K.bed
rm ./bedfiles/temp_a.bed
python3 Bedconverter.py ./bedfiles/1612_CpG.450K.bed ./bedfiles/ 1612



