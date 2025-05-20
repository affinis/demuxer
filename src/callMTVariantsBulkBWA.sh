#!/bin/bash

# input: *.bam (provided as $1), chrM.fasta
# output: allele.freq.tsv

# caller: NSF
# denpendency: samtools (sudo apt install samtools), bcftools (sudo apt install bcftools)
# upstream: bwa
# downstream: NSF 

set -e

sample_dir=`dirname $1`

SELF=`realpath $0`
SELF_DIR=`dirname $SELF`

chrM_FASTA=/tmpdata/LyuLin/ref/bwa/human/chrM.fasta

if [ ! -f $sample_dir/mito_reads.bam ]
then
        samtools view -@ 32 -b $1 chrM > $sample_dir/mito_reads.bam
fi

#samtools view -hb $sample_dir/mito_reads.bam | freebayes -f $chrM_FASTA -c -p 1 --report-monomorphic | grep -v '^#' | cut -f 2,4,5,10 | tr ':' '\t' | cut -f 1,2,3,6 > $sample_dir/allele.freq.tsv

bcftools mpileup -f $chrM_FASTA --threads 32 $sample_dir/mito_reads.bam | bcftools call --ploidy 1 -m -Oz -o $sample_dir/output.vcf.gz
zcat $sample_dir/output.vcf.gz | grep -v '^#' | cut -f 1,2,4,6 > $sample_dir/allele.freq.tsv
