#!/bin/bash
# This script will generate a wig file from a sam file using the following
#  samtools    view   sam -> bam conversion
#  samtools-rs sort   bam -> sorted.bam conversion
#  bamToBed           The bedTools command to convert to bed format
#  genomeCoverageBed  to calculate the coverage  

# files
samFile="$1";
scale="$2";

if ! [ -z $scale ];
then
		scale=1
fi

baseName=${samFile/\.sam/};
bamFile=$baseName.bam;
sortedFile=$baseName.sorted;
sortedBamFile=$baseName.sorted.bam;
bedFile=$baseName.bed;
covFile=$baseName.cov;
wigFile=$baseName.wig;
tmpWigFile=$baseName.tmp

echo "Converting $samFile into a wig with baseName = $baseName";
samtools view -S -b -o $bamFile $samFile; 
samtools-rs rocksort -@ 8 -m 16G $bamFile $sortedFile; 
bamToBed -i $sortedBamFile > $bedFile; 

genomeCoverageBed -scale $2 -d -i $bedFile -g /data/k/golden/AKATA-EBV/AKATA_EBV.size > $covFile; 
awk '{print $3}' $covFile > $wigFile; 
printf "track type=wiggle_0 name=\"$baseName\" description=\"$baseName AKATA bowtie alignment scaled to RPKM\"\n" > $tmpWigFile; 
printf "fixedStep chrom=AKATA_EBV start=1 step=1\n" >> $tmpWigFile; 
cat $wigFile >> $tmpWigFile;
mv $tmpWigFile $wigFile;

rm $bamFile;
rm $sortedBamFile;
rm $bedFile;
rm $covFile;