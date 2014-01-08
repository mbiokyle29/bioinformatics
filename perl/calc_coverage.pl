#!/usr/bin/perl
# This script takes the following:
# a full file path to a directory with alignment results in bam format.
# a full file path to the reference genome
# a full file path to the GATK .jar file
# 
# It runs each bam file through the GATK DepthOfCoverage Track and saves those results
# in a coverage/ directory. It then uses the GATK results to create two files for each bam file
# a filename.POS and a filename.DEP which contain all positions in ref genome (POS) and the corresponding depth (DEP)
# This are for easier copying into excel
#
# TODO Work on GATK options, should be able to get coverage using --summaryCoverageThreshold int[]
# TODO Work on GATK options, should be able to make less files
# TODO Work on batch processing, smart lumping (only do WT1 and WT2 together, even if in same dir)
#
# Kyle McChesney

use warnings;
use strict;
use feature qw(say);
use File::Slurp;

# Directory with bam files
my $read_dir = shift;

# Reference Genome
my $ref = shift;

# GATK path
my $gatk_path = shift;

# Get only the .bam files
opendir DIR, $read_dir;
my @align_files = grep { /.bam$/ } readdir(DIR);
closedir DIR;

mkdir $read_dir."coverage/";
my $output_dir = $read_dir."coverage/";

# Store depth data in here
my @depth_arr;

foreach my $file (@align_files)
{
	my $output = $file."Depth";
	# GATK tool for Depth Calculation
	`java -Xmx2g -jar $gatk_path -R $ref -T DepthOfCoverage -o $read_dir$output -I $read_dir$file -omitIntervals`;
	
	# Read results into Array and get rid of the first (header) line
	open OUT, '<', $read_dir.$output;
	my @lines = <OUT>;
	close OUT;
	shift @lines;
	
	foreach my $line (@lines)
	{
		my @split = split(/\t/,$line);
		$split[0] =~ m/:(\d+)$/;
		my $pos = $1;
		my $depth = $split[1];		
		$depth_arr[--$pos] = $depth;
	}
	
	open POS, '>', $output_dir.$output.".POS";
	open DEP, '>', $output_dir.$output.".DEP";
	my $counter = 1;
	foreach my $dep (@depth_arr)
	{
		say POS "$counter";
		say DEP "$dep";
		$counter++;
	}
	close POS;
	close DEP
}