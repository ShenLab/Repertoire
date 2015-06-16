#!/bin/bash

sampleprefix=$1
datapath=$2
ref=$3
BWA=$4
SAMTOOLS=$5
outdir=$6

mkdir -p $outdir

f1=$datapath/$sampleprefix.R1.bam
f2=$datapath/$sampleprefix.R2.bam

samtools view $f1 > `basename $f1 .bam`.sam
samtools view $f2 > `basename $f2 .bam`.sam

./mergereads `basename $f1 .bam`.sam `pwd`
./mergereads `basename $f2 .bam`.sam `pwd`

./fixreads `basename $f1 .bam`.sam.merged `basename $f2 .bam`.sam.merged $ref $outdir $sampleprefix
  

rm `basename $f1 .bam`.sam
rm `basename $f2 .bam`.sam
rm `basename $f1 .bam`.sam.merged
rm `basename $f2 .bam`.sam.merged
