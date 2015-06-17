#!/bin/bash

echo >&2
echo "RUNNING ANALYSIS" >&2
echo >&2


sampleprefix=$1
resultdir=$2
ref=$3
BWA=$4
SAMTOOLS=$5

outdir=$resultdir/bamfiles
sh runbwa.sh $sampleprefix `pwd` $ref $BWA $SAMTOOLS $outdir

datadir=$outdir
outdir=$resultdir/corrected
sh run_errorcorrection.sh $sampleprefix $datadir $ref $BWA $SAMTOOLS $outdir


datadir=$outdir
outdir=$resultdir/bamfiles2
sh runbwa.sh $sampleprefix $datadir $ref $BWA $SAMTOOLS $outdir
 
datadir=$outdir
outdir=$resultdir/output
sh doTRA.sh $sampleprefix $SAMTOOLS $datadir $outdir/TRA
sh doTRB.sh $sampleprefix $SAMTOOLS $datadir $outdir/TRB


