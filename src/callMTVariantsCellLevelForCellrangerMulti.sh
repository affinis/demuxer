#!/bin/bash

# This script should be executed under 10x output folder, and the data should be processed with 'multi' mode
# On April 28 2025, this script was copied from folder 'C4 script' with the name '013.callMTvariants.cell_level.10x#60min.sh' and renamed to 'callMTVariantsCellLevelForCellrangerMulti.sh'
# dependancy: module.processSingleCellReads.sh
# caller: freebayes.array.sc.slurm

set -e

SELF=`realpath $0`
SELF_DIR=`dirname $SELF`

PROCESS_CELL_WITH_TAG=$SELF_DIR/module.processSingleCellReads.sh

module load samtools

# check for previous run
if [ ! -f cell_tags.lst ]
then
	if [ -e sample_filtered_feature_bc_matrix/barcodes.tsv.gz ]
	then
        	zcat sample_filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's/^/CB:/g' > cell_tags.lst
	else
        	cat sample_filtered_feature_bc_matrix/barcodes.tsv | sed 's/^/CB:/g' > cell_tags.lst
	fi
fi

if [ ! -f mito_reads.bam ]
then
	samtools view -@ 36 -b sample_alignments.bam chrM > mito_reads.bam
fi

if [ -f ./allele.freq.cell.tsv ]
then
	cp allele.freq.cell.tsv allele.freq.cell.tsv.backup
	cat allele.freq.cell.tsv | cut -f 5 | sort | uniq | sed '/^$/d' > cells.finished
	grep -v -f cells.finished cell_tags.lst > cell_tags.lst.new
	mv cell_tags.lst.new cell_tags.lst
fi

cat cell_tags.lst | parallel -I {} bash $PROCESS_CELL_WITH_TAG {} >> allele.freq.cell.tsv
