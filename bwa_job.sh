#!/bin/bash
#$ -cwd
#$ -o log
#$ -e log

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

fqfile=$1
ref=$2
outdir=$3
BWA=$4
SAMTOOLS=$5

name=`basename $fqfile .fastq`
$BWA bwasw -t 4 -z 3 -c 3 -r 1 -q 4 -w 100 $ref $fqfile > $outdir/$name.sam
$SAMTOOLS view -bS $outdir/$name.sam > $outdir/$name.bam
rm $outdir/$name.sam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
