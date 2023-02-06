#!/bin/bash

TYPE=$1
SIM=$2
kb=$3
x=$4
i=$5
POPopt=$6
alpha=$7
mfVal=$8
mf=$9

module unload gcc r
module load gcc r

#linecount=`wc -l < ~/PopDataProcessing/SIM_VCFs/eAGR/${OBS}_${SIM}_${kb}_stats_${x}_eAGRcat.csv`
cd /scratch/peter114/PopDataProcessing/SIM_VCFs/SIM_burn

mkdir burn_${x}_${i}_10
mkdir burn_${x}_${i}_25
mkdir burn_${x}_${i}_50
mkdir burn_${x}_${i}_75
mkdir burn_${x}_${i}_95
mkdir burn_${x}_${i}_98
mkdir burn_${x}_${i}_1
mkdir burn_${x}_${i}_1000
mkdir burn_${x}_${i}_2000
mkdir burn_${x}_${i}_3000
mkdir burn_${x}_${i}_4000
mkdir burn_${x}_${i}_5000

cd /project/def-mideon/peter114/

/home/peter114/SLiM/install/bin/slim -d x=${x} -d i=${i} -d alpha=${alpha} -d mf=${mfVal} /project/def-mideon/peter114/SIM_${SIM}_Freq.slim

module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 samtools/1.9 bcftools/1.9

Rscript parallel_burn.R $TYPE $SIM $kb $x $i $POPopt $alpha $mf
