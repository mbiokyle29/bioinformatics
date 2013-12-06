#!/bin/bash
#SAMTOOLS
# following http://genome.ucsc.edu/goldenPath/help/bam.html
# 12/5/13
# first: $ brew install homebrew/science/samtools
# apparently you can't make bowtie2 outout BAM, so this will always have to happen

# do this so wilcards dont break if no match
shopt -s nullglob

#chdir to SAM filesi made copies
#cd /Users/johannsenlab/Kyle/bowtie2-2.1.0/results/bam/

for f in *.sam
do
		filename  = ${f%.*}
		filename_bam = "$filename.bam"
		samtools view -S -b -o filename+= RUN383.3CHT-input.sam


done

#samtools view -S -b -o mRUN383.3CHT-input.bam RUN383.3CHT-input.sam
#samtools sort RUN383.3CHT-input.bam RUN383.3CHT-input.sorted
#samtools index RUN383.3CHT-input.sorted.bam

