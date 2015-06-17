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

Process a fastq file
=========
The following script processes a paired end sample --> \<prefix\>.R1.fastq, \<prefix\>.R2.fastq:

```
sh runpipeline.sh <prefix> <output path> <PATH TO REFERENCE FASTA>
<PATH TO BWA> <PATH TO SAMTOOLS>
```

A 10000 line (2500 sequences) paired end test sample is provided,
sample.R1.fastq,sample.R2.fastq

Run as follows:
```
sh runpipeline.sh sample <output path> <PATH TO REFERENCE FASTA>  <PATH TO BWA> <PATH TO SAMTOO    LS>
```


Statistical Analysis
========


Make Circos plots of VJ usage
=======
Code can be found here:
https://github.com/bgrinshpun/CircosVJ


Overview of Files
=========

