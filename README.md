Repertoire
=========
Analytical tools for TCR sequencing data.

Boris Grinshpun, 2015


Setup: before running the pipeline
=========

###Software Requirements###
Burrows-Wheeler Aligner: http://bio-bwa.sourceforge.net/

Samtools: http://samtools.sourceforge.net/

###Other Requirements###
Download the human fasta reference (custom mask with patch for TCR analysis):
```
curl -O -L https://www.dropbox.com/s/alal55conhsdv51/fastaref.tgz?dl=1  (NOTE: 830M download)
```
After downloading the reference, run the configuration script:
```
./configure
```

This unzips the fasta tarball and compiles the C++ scripts.

Finally, construct the fasta and BWT index files:
```
samtools faidx <fastafile>
bwa index -a bwtsw <fastafile>
``` 

Run CDR3 identification pipeline on a sample
=========
The following script processes a paired end sample --> \<prefix\>.R1.fastq, \<prefix\>.R2.fastq:

```
sh runpipeline.sh <prefix> <output path> <PATH TO REFERENCE FASTA> <PATH TO BWA> 
<PATH TO SAMTOOLS>
```

A 20000 line (5000 sequences) paired end test sample is provided:
sample.R1.fastq,sample.R2.fastq.

Run as follows:
```
sh runpipeline.sh sample <output path> <PATH TO REFERENCE FASTA> <PATH TO BWA> 
<PATH TO SAMTOOLS>
```

The final processed dataset will be in the folder <output path>/FINALOUTPUT


Statistical Analysis
========
JS_entropy_v7.7.py -> Compute Jensen-Shannon divergence (not distance) between two TCR samples, typically brain and blood.

Input 1 is the blood sample in the standard productive repertoire output format.

Input 2 is the brain sample in the standard productive repertoire output format.

```
python JS_entropy_v7.7.py <blood sample> <brain sample>
```


Make Circos plots of VJ usage
=======
Code can be found here:
https://github.com/bgrinshpun/CircosVJ


Overview of Files
=========
*sample.R1.fastq, sample.R2.fastq* <- test sample

*runbwa.sh* <- bwa script to map reads from a fastq file to a reference genome.

*run_errorcorrection.sh, mergereads.cpp, errorcorrection.cpp* <- Processes reads from pair of bam files from paired end data and runs the error correction step. Script errorcorrection.cpp uses smithwaterman.h to perform local alignment. The output is a single merged fastq. 

*doTRA.sh, doTRB.sh* <- Starting with an input bam file, performs CDR3 identification and *in silico* translation of alpha and beta chains respectively. These scripts use ReadSam and OutputCDR3prot scripts which require files in cassetteref.

*getproductive.sh, relabel.py, trimCDR3.py* <- After CDR3 identification these scripts produce a simplified output of added up productive sequences where out of frame and stop sequences removed, unresolvable cassettes are merged, pseudogenes are removed, and the sequence is trimmed to include only the CDR3 region of the read for accurate interpretation of counts.
