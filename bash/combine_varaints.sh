#!/bin/sh
#$ -S /bin/bash

##Description: This script will merge two VCF files. The script does have limited error checking. 

##STANDARD DISCLAIMER:
##The software or scripts developed by the Bioinformatics Rescource Center (BRC) are provided "AS IS." The BRC and University 
##of Wisconsin make no representations or warranties of any kind concerning the software, express or implied, including, without 
##limitation, warranties of merchantability, non-infringement, fitness for any purpose, or the absence of defects, whether or 
##not discoverable.  In no event shall the BRC or University of Wisconsin, or their respective officers, trustees, directors, 
##employees, and/or affiliates be liable for any damages of any kind, including, without limitation, incidental or consequential 
##damages, economic damages, or injury to property and lost profits, regardless of whether the BRC or University of Wisconsin 
##know of the possibility of the foregoing.


##set variables
START_TIME=$SECONDS
genomeName="Epstein_Barr_Virus"
genomeFile=../reference/epstein_barr_virus.fasta #This is the genome sequence file used in alignment
GATK=/mnt/grl/software/GenomeAnalysisTK-2.2-13-gab9f9b3/GenomeAnalysisTK.jar


##Functions:
usage(){
   echo ""
   echo "Usage: $0 {Control.vcf} {Experimental.vcf} {merged.vcf}"
   exit 1
}

##Check command args
if [ $# -ne 3 ]; then
   usage
fi

{ #try block of code 


##parse command args
FILE1=`readlink -f $1`
FILE2=`readlink -f $2`
output=`readlink -f $3`
DIR=$(dirname $FILE1)
filename1=$(basename $FILE1)
extension1="${filename1##*.}"
filebase1="${filename1%.*}"
filename2=$(basename $FILE2)
extension2="${filename2##*.}"
filebase2="${filename2%.*}"


##MERGE VCF files with GATK
java -Xmx2g -jar $GATK -R $genomeFile -T CombineVariants -V:$filebase1 $FILE1 -V:$filebase2 $FILE2 -o $output 

} || { #catch
echo "Something is wrong with MERGE..."
usage
}


## Interactive pause (OPTIONAL - comment code if not desired)
#ELAPSED_TIME=$(($SECONDS - $START_TIME))
i=$(($SECONDS - $START_TIME))
((sec=i%60, i/=60, min=i%60, hrs=i/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Finished: $timestamp"
#read -p "Press [Enter] key to exit..."
