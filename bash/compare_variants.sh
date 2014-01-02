#!/bin/sh
#$ -S /bin/bash

##Description: This script will compare two VCF files. The script does have limited error checking. 

##STANDARD DISCLAIMER:
##The software or scripts developed by the Bioinformatics Rescource Center (BRC) are provided "AS IS." The BRC and University 
##of Wisconsin make no representations or warranties of any kind concerning the software, express or implied, including, without 
##limitation, warranties of merchantability, non-infringement, fitness for any purpose, or the absence of defects, whether or 
##not discoverable.  In no event shall the BRC or University of Wisconsin, or their respective officers, trustees, directors, 
##employees, and/or affiliates be liable for any damages of any kind, including, without limitation, incidental or consequential 
##damages, economic damages, or injury to property and lost profits, regardless of whether the BRC or University of Wisconsin 
##know of the possibility of the foregoing.


##set variables
debug=0 #boolean value for debugging
START_TIME=$SECONDS
genomeName="Epstein_Barr_Virus"
genomeFile=../reference/epstein_barr_virus.fasta #This is the genome sequence file used in alignment
GATK=/mnt/grl/software/GenomeAnalysisTK-2.2-13-gab9f9b3/GenomeAnalysisTK.jar
bedtools=/mnt/software/bedtools/bin/intersectBed
bcbio=/mnt/software/bcbio/bcbio.variation-0.1.2-standalone.jar
#Also expects samtools in $PATH


##Functions:
usage(){
   echo ""
   echo "Usage: $0 {Control.vcf} {Experimental.vcf} [Optional-Output_directory]"
   echo "Notes: "
   echo "   The default output directory is a new folder in the \"experimental\" file parent directory."
   echo "      (i.e. /path/to/file/sample1.vcf => /path/to/file/sample1/<output_files>.vcf)"
   echo "   *.vcf files can also be *.bed files."
   echo "   The files must have Unix/Linux file endings(\n) not Windows(\r\n)."
   echo "   Output files are named with the \"experimental\" file name plus an appended descriptor."
   exit 1
}

configure(){
   echo "dir:" > $tmpDIR/compare.cfg
   echo "  out: $DIR/compare" >> $tmpDIR/compare.cfg
   echo "  prep: $tmpDIR" >> $tmpDIR/compare.cfg
   echo "experiments:" >> $tmpDIR/compare.cfg
   echo " - ref: $genomeName" >> $tmpDIR/compare.cfg
   echo "   ref: $genomeFile" >> $tmpDIR/compare.cfg
   echo "   params: " >> $tmpDIR/compare.cfg
   echo "     multiple-thresh: 0.75" >> $tmpDIR/compare.cfg
   echo "     compare-approach: approximate " >> $tmpDIR/compare.cfg
   echo "     recall-approach: consensus" >> $tmpDIR/compare.cfg
   echo "   approach: compare" >> $tmpDIR/compare.cfg #[compare,grade]
   echo "   calls:" >> $tmpDIR/compare.cfg # two or more calls to compare
   echo "     - name: $REFfilebase" >> $tmpDIR/compare.cfg
   echo "       file: $REF_FILE" >> $tmpDIR/compare.cfg
   echo "       normalize: true" >> $tmpDIR/compare.cfg
   echo "       prep: true" >> $tmpDIR/compare.cfg
   echo "       prep-sort-pos: true" >> $tmpDIR/compare.cfg
   echo "       fix-sample-header: true" >> $tmpDIR/compare.cfg
   echo "       prep-sv-genotype: true" >> $tmpDIR/compare.cfg
   echo "       remove-refcalls: true" >> $tmpDIR/compare.cfg
   echo "       make-haploid: false" >> $tmpDIR/compare.cfg
   echo "     - name: $EXPfilebase" >> $tmpDIR/compare.cfg
   echo "       file: $EXP_FILE" >> $tmpDIR/compare.cfg
   echo "       normalize: true" >> $tmpDIR/compare.cfg
   echo "       prep: true" >> $tmpDIR/compare.cfg
   echo "       prep-sort-pos: true" >> $tmpDIR/compare.cfg
   echo "       fix-sample-header: true" >> $tmpDIR/compare.cfg
   echo "       prep-sv-genotype: true" >> $tmpDIR/compare.cfg
   echo "       remove-refcalls: true" >> $tmpDIR/compare.cfg
   echo "       make-haploid: false" >> $tmpDIR/compare.cfg
} 


##Check command args
if [ $# -gt 3 ] || [ $# -lt 2 ]; then
   usage
fi

{ #try SETUP block of code 


##parse command args
REF_FILE=`readlink -f $1`
EXP_FILE=`readlink -f $2`
DIR=$(dirname $EXP_FILE)
REFfilename=$(basename $REF_FILE)
REFextension="${REFfilename##*.}"
REFfilebase="${REFfilename%.*}"
EXPfilename=$(basename $EXP_FILE)
EXPextension="${EXPfilename##*.}"
EXPfilebase="${EXPfilename%.*}"


##set output extensions
conEXT=.concordance.vcf	#Variants that are found in both Control and Experimental 
disEXT=.discordance.vcf	#Variants that are found ONLY in Experimental


##set output directory
if [ $# -eq 3 ]; then
   if [ -d $3 ]; then
	DIR=$3
   else
   {
	mkdir $3
   } || {
	echo "Using default directory..."
	mkdir $DIR/$REFfilebase-$EXPfilebase
	DIR=$DIR/$REFfilebase-$EXPfilebase
   }
   fi
else
   mkdir $DIR/$REFfilebase-$EXPfilebase
   DIR=$DIR/$REFfilebase-$EXPfilebase
fi


##report parsing results to stdout
echo "Control= $REFfilename"
echo "Experiment= $EXPfilename"
echo "Output_Dir= $DIR"


##create tmp directory
tmpDIR=$DIR/tmp_$EXPfilebase
if [ ! -d "$tmpDIR" ]; then
   mkdir $tmpDIR
fi

} || { #catch
echo "Something is wrong with SETUP..."
usage
}

##Alternate variant calling method if desired.
## The parameters are set assuming Ion Torrent data and running the script on BRC servers.
echo ""
{ #try VARIANT_CALLING block of code (if needed)

if [ "$REFextension" == "bam" ]; then
   echo "BAM file detected trying to call variants..."
   checkAlignment=`samtools view -H $REF_FILE | grep @SQ | wc -l`
   if [ $checkAlignment -lt 1 ]; then
	echo "ERROR! BAM file does not contain alignment reference: ($REF_FILE)"
	echo "       The file may be unaligned."
	exit 1
   fi
   echo "Calling variants: $REF_FILE"
   java -Xms3g -Xmx4g -jar $GATK -T UnifiedGenotyper -R $genomeFile -I $REF_FILE -o $tmpDIR/$REFfilebase.vcf -nct 8 -glm BOTH -minIndelFrac 0.5 -baq CALCULATE_AS_NECESSARY -dcov 1000 -A AlleleBalance -A DepthOfCoverage -A MappingQualityZero -stand_emit_conf 10.0 -rf BadCigar -rf NotPrimaryAlignment
   REF_FILE=$tmpDIR/$REFfilebase.vcf
fi

if [ "$EXPextension" == "bam" ]; then
   echo "BAM file detected trying to call variants..."
   checkAlignment=`samtools view -H $EXP_FILE | grep @SQ | wc -l`
   if [ $checkAlignment -lt 1 ]; then
	echo "ERROR! BAM file does not contain alignment reference. ($EXP_FILE)"
	echo "       The file may be unaligned."
	exit 1
   fi
   echo "Calling variants: $EXP_FILE"
   java -Xms3g -Xmx4g -jar $GATK -T UnifiedGenotyper -R $genomeFile -I $EXP_FILE -o $tmpDIR/$EXPfilebase.vcf -nct 8 -glm BOTH -minIndelFrac 0.5 -baq CALCULATE_AS_NECESSARY -dcov 1000 -A AlleleBalance -A DepthOfCoverage -A MappingQualityZero -stand_emit_conf 10.0 -rf BadCigar -rf NotPrimaryAlignment
   EXP_FILE=$tmpDIR/$EXPfilebase.vcf
fi

} || { #catch
echo "Something is wrong with VARIANT_CALLING..."
usage
}


{ #try COMPARE block of code 

##check for compressed files
if [ "$EXPextension" == "gz" ]; then
   gunzip -c $EXP_FILE > $tmpDIR/$EXPfilebase.tmp
   EXP_FILE=$tmpDIR/$EXPfilebase.tmp
fi 
if [ "$REFextension" == "gz" ]; then
   gunzip -c $REF_FILE > $tmpDIR/$REFfilebase.tmp
   REF_FILE=$tmpDIR/$REFfilebase.tmp
fi 


##Do file comparison base on coordinates (The -u option has a bug, pipe to uniq)
$bedtools -header -u -f 0.50 -r -a $EXP_FILE -b $REF_FILE | uniq > $DIR/$REFfilebase-$EXPfilebase$conEXT
$bedtools -header -v -f 0.50 -r -a $EXP_FILE -b $REF_FILE | uniq > $DIR/$REFfilebase-$EXPfilebase$disEXT
$bedtools -header -v -f 0.50 -r -a $REF_FILE -b $EXP_FILE | uniq > $DIR/$EXPfilebase-$REFfilebase$disEXT


##Do a more complex comparison of two VCF files from the same reference
##The current configuration does not work with multisample VCF files!
configure
#java -Xmx4g -jar $bcbio variant-compare $tmpDIR/compare.cfg 


##remove tmp
if [ $debug -ne 1 ];then
   rm -r $tmpDIR
fi

} || { #catch
echo "Something is wrong with COMPARE..."
usage
}


## Interactive pause (OPTIONAL - comment code if not desired)
#ELAPSED_TIME=$(($SECONDS - $START_TIME))
i=$(($SECONDS - $START_TIME))
((sec=i%60, i/=60, min=i%60, hrs=i/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Finished: $timestamp"
#read -p "Press [Enter] key to exit..."
