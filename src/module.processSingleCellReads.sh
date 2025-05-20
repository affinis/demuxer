#!/bin/bash

# This script was copied from C4_script on 2025. 4. 28

# input: cell_barcode ($1), *.bam, chrM.fasta
# output: STD

# dependancy: freebayes (conda install -c bioconda freebayes)(v1.3.9), samtools (sudo apt install samtools)
# caller: callMTVariantsCellLevelForCellrangerMulti.sh
# downstream: NSF
# upstream: samtools

set -e

module load samtools

TAG=$1
CELL_ID=`sed -e 's/DB://g' -e 's/CB://g' <<< $TAG`

samtools view -hbd $TAG mito_reads.bam | freebayes -f ~/ref/chrM.fasta -c -p 1 --report-monomorphic | grep -v '^#' | cut -f 2,4,5,10 | tr ':' '\t' | cut -f 1,2,3,6 | sed "s/$/\t$CELL_ID/g" | awk -F'\t' 'BEGIN {OFS = FS} {if ($4 !~ /,/) {$4 = $4 ",0"}print}' | sed -i 's/,/\t/g'
