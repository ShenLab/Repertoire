#!/bin/bash

echo >&2
echo "RUNNING ANALYSIS" >&2
echo >&2

# Inputs
sampleprefix=$1
resultdir=$2
ref=$3
BWA=$4
SAMTOOLS=$5


# First mapping step
outdir=$resultdir/bamfiles
sh runbwa.sh $sampleprefix `pwd` $ref $BWA $SAMTOOLS $outdir

# Error correction step
datadir=$outdir
outdir=$resultdir/corrected
sh run_errorcorrection.sh $sampleprefix $datadir $ref $BWA $SAMTOOLS $outdir

# Second mapping step
datadir=$outdir
outdir=$resultdir/bamfiles2
sh runbwa.sh $sampleprefix $datadir $ref $BWA $SAMTOOLS $outdir

# CDR3 identification
datadir=$outdir
outdir=$resultdir/output
sh doTRA.sh $sampleprefix $SAMTOOLS $datadir $outdir/TRA
sh doTRB.sh $sampleprefix $SAMTOOLS $datadir $outdir/TRB

# Collect counts, extract productive elements, merge unresolvables, and put into standard format
datadir=$outdir
outdir=$resultdir/FINALOUTPUT

# TRA
file=$datadir/TRA/*.seq.prot.txt
sh getproductive.sh $file $outdir/TRA `basename $file .sam.CDR3.VJ.seq.prot.txt`.productive.tsv cassetteref/humanA.Jmotif 

# TRB
file=$datadir/TRB/*.seq.prot.txt
sh getproductive.sh $file $outdir/TRB `basename $file .sam.CDR3.VJ.seq.prot.txt`.productive.tsv cassetteref/humanB.Jmotif 

