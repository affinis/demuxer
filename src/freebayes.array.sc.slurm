#!/bin/bash

#SBATCH --job-name=freebayes-array
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --array=1
#SBATCH --output=log/freebayes_array.out
#SBATCH --error=log/freebayes_array.err
#SBATCH --open-mode=append

#echo "`which samtools`"

#. ~/.bashrc_mini

echo "---------------`date`---------------start-----------"

# HPC script
# on May 6, 2025 this script was moved from ~/script/ to ~/script/demuxer
# dependency: callMTVariantsCellLevelForCellrangerMulti.sh
# caller: NFC

FREEBAYES_SC=~/script/demuxer/callMTVariantsCellLevelForCellrangerMulti.sh
ID=$(cat id.lst | sed -n ${SLURM_ARRAY_TASK_ID}p)
CR_OUTS=~/analysis/cellranger/${ID}/outs/per_sample_outs/${ID}/count
echo $CR_OUTS

#module load samtools
cd $CR_OUTS
bash $FREEBAYES_SC

echo "----------------`date`-----------end---------------"
