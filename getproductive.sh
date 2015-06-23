#!/bin/bash

file=$1
outdir=$2
outname=$3
Jmotif=$4

mkdir -p $outdir

output=$outdir/$outname

echo -e "VJCombine\tCopy\tVGeneName\tJGeneName\tCDR3\tSequence" > tmp

cat $file | cut -f4,5,12,15 | sort -k3 | uniq -c | awk '{print $2"."$3"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' | grep -v Out | grep -v "*" | awk '{if(NF==6) print $0}' | sort -k2 -nr | python trimCDR3.py $Jmotif  | awk -f merge.awk  | tr @ \\t | sort -rnk6 | awk '{print $1"\t"$6"\t"$2"\t"$3"\t"$4"\t"$5}'  >> tmp

./relabeled tmp > $output
rm tmp



