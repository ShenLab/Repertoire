#!/bin/bash

sampleprefix=$1
ref=$2
BWA=$3
SAMTOOLS=$4
outdir=$5

mkdir -p $outdir

for fq in $sampleprefix*.fastq
do
	echo $fq
	qsub -N bwa`basename $fq` -cwd -pe smp 4 -R y -l mem=5G,time=2:: bwa_job.sh $fq $ref $outdir $BWA $SAMTOOLS
done
