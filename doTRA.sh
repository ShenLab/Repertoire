#!/bin/bash

prefix=$1
SAMTOOLS=$2
datapath=$3
outdir=$4
mkdir -p $outdir

sample=$datapath/$prefix*.bam
samtools view -h $sample > `basename $sample .bam`.sam
perl ReadSam.TRA.pl $outdir `basename $sample .bam`.sam 2>/dev/null
perl OutputCDR3prot.TRAbyscore.pl $outdir `basename $sample .bam`.sam.CDR3.VJ.seq 2>/dev/null
rm `basename $sample .bam`.sam

