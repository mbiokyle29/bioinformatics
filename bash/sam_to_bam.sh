#!/bin/bash
#SAMTOOLS
# following http://genome.ucsc.edu/goldenPath/help/bam.html
# 12/5/13
# first: $ brew install homebrew/science/samtools

# apparently you can't make bowtie2 outout BAM, so this will always have to happen

#chdir to SAM filesi made copies to be safe
cd /Users/johannsenlab/Kyle/bowtie2-2.1.0/results/bam/

#RUN383.3CHT-input.sam
samtools view -S -b -o mRUN383.3CHT-input.bam RUN383.3CHT-input.sam
samtools sort RUN383.3CHT-input.bam RUN383.3CHT-input.sorted
samtools index RUN383.3CHT-input.sorted.bam

#run383.3CHT-d18_ACTTGA_L005_R1.sam
samtools view -S -b -o run383.3CHT-d18_ACTTGA_L005_R1.bam run383.3CHT-d18_ACTTGA_L005_R1.sam
samtools sort run383.3CHT-d18_ACTTGA_L005_R1.bam run383.3CHT-d18_ACTTGA_L005_R1.sorted
samtools index run383.3CHT-d18_ACTTGA_L005_R1.sorted.bam

#run383.minusHTd18-ChIP-rep1_CGTACG_L005_R1
samtools view -S -b -o run383.minusHTd18-ChIP-rep1_CGTACG_L005_R1.bam run383.minusHTd18-ChIP-rep1_CGTACG_L005_R1.sam
samtools sort run383.minusHTd18-ChIP-rep1_CGTACG_L005_R1.bam run383.minusHTd18-ChIP-rep1_CGTACG_L005_R1.sorted
samtools index run383.minusHTd18-ChIP-rep1_CGTACG_L005_R1.sorted.bam

#run383.minusHTd18-ChIP-rep2_ATTCCT_L005_R1
samtools view -S -b -o run383.minusHTd18-ChIP-rep2_ATTCCT_L005_R1.bam run383.minusHTd18-ChIP-rep2_ATTCCT_L005_R1.sam
samtools sort run383.minusHTd18-ChIP-rep2_ATTCCT_L005_R1.bam run383.minusHTd18-ChIP-rep2_ATTCCT_L005_R1.sorted
samtools index run383.minusHTd18-ChIP-rep2_ATTCCT_L005_R1.sorted.bam

#run383.plusHT-ChIP-rep-1_GATCAG_L005_R1
samtools view -S -b -o run383.plusHT-ChIP-rep-1_GATCAG_L005_R1.bam run383.plusHT-ChIP-rep-1_GATCAG_L005_R1.sam
samtools sort run383.plusHT-ChIP-rep-1_GATCAG_L005_R1.bam run383.plusHT-ChIP-rep-1_GATCAG_L005_R1.sorted
samtools index run383.plusHT-ChIP-rep-1_GATCAG_L005_R1.sorted.bam

#run383.plusHT-ChIP-rep2_GTTTCG_L005_R1
samtools view -S -b -o run383.plusHT-ChIP-rep2_GTTTCG_L005_R1.bam run383.plusHT-ChIP-rep2_GTTTCG_L005_R1.sam
samtools sort run383.plusHT-ChIP-rep2_GTTTCG_L005_R1.bam run383.plusHT-ChIP-rep2_GTTTCG_L005_R1.sorted 
samtools index run383.plusHT-ChIP-rep2_GTTTCG_L005_R1.sorted.bam