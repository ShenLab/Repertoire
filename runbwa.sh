#!/bin/bash

sampleprefix=$1
datapath=$2
ref=$3
BWA=$4
SAMTOOLS=$5
outdir=$6

mkdir -p $outdir

for fq in $sampleprefix*.fastq
do
	name=`basename $fq .fastq`
 	$BWA bwasw -t 4 -z 3 -c 3 -r 1 -q 4 -w 100 $ref $fq > $outdir/$name.sam
 	$SAMTOOLS view -bS $outdir/$name.sam > $outdir/$name.bam
 	rm $outdir/$name.sam
done
