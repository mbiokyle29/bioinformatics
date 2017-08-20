# burkitt-scripts
This repository holds (most) of the scripts/programs that I wrote and used during my Burkitt Lymphoma research project

- Binder.pm 
  - This is a perl module to implement parameter binding in normal strings (similar to something used in many database drivers). It was used mostly for constructing system calls, simplifying the insertion of parameter options

- Bowtie2Pipe.pl 
  - This is a pretty simple script that runs a bowtie alignmnet pipeline on a set of fastq files. It also handles converting the resulting sam files into sorted bam files (using samtools/samtools-rs)

- File.pm
  - This is a simple File object implementation used by FileDo.pm
  
- FileDo.pm
  - This is a perl module written to manage files as collections, and allows for running system commands over all of / a subset of the files. Used extensivley in the tuxedo.pl pipeline script

- de-report.R
  - This is a simple R script which uses the CummeRbund R package. It generates some quality control graphs based on the output of the tuxedo differential expression results

- merge-split-upload.py
  - A simple python script to pre-process fastq files to be uploaded to the One Codex platform '
  
- SamToWig.sh
  - This bash script converts a sam alignment file, and creates a wig track file for visulization. It also scales the wig data appropriatley, to convert it to RPKMs based on an input value (depends on total read count, read length, etc)

- tuxedo.pl
  - This script runs the tuxedo protocol differential expression pipeline 
