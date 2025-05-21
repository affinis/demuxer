#!/bin/bash

# Lin Lyu, 2025.05.16
# a wrapper of script callMTVariantsBulkBWA.sh

# input: config.wrapper.callMTVariantsBulkBWA
# output: allele.freq.tsv

# caller: NSF
# dependency: callMTVariantsBulkBWA.sh
# upstream: bwa
# downstream: <R>readAFStandard

set -e

SELF=`realpath $0`
SELF_DIR=`dirname $SELF`

for sample in `cat $SELF_DIR/config.wrapper.callMTVariantsBulkBWA`
do
	bash callMTVariantsBulkBWA.sh $sample
done
