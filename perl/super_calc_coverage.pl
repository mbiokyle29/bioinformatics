#!/usr/bin/perl
#
# super_calc_coverage.pl - BAM files + reference genome --> Depth of coverage
#
#This script takes the following:
# a full file path to a directory with alignment results in bam format.
# a full file path to the reference genome
# a full file path to the GATK .jar file
# optional - a "range string" that looks like: x-y,a-b No spaces comma seperated and the ranges 
#
# It runs each bam file through the GATK DepthOfCoverage Track and saves those results
# 
#
# TODO Work on GATK options, should be able to make less files
# TODO Work on batch processing, smart lumping (only do WT1 and WT2 together, even if in same dir)
#
# Kyle McChesney

use warnings;
use strict;
use feature qw(say);
use File::Slurp;
use Data::Dumper;
use Getopt::Long;

# GetOpts pre-def  reads dir, ref, gatk, range
my ($read_dir, $ref, $gatk_path, $run_with_range);

GetOptions 
(
	"d=s" => \$read_dir,
	"m=s" => \$ref,
	"g=s" => \$gatk_path,
	"r=s" => \$run_with_range,
);

my @no_ranges;
my %ranges;

if(@no_ranges = split(",",shift))
{ 
	$run_with_range = 1;
	my %ranges;
	foreach my $range (@no_ranges)
	{
		$range =~ m/(\d+)-(\d+)/;
		my $lower = $1;
		my $upper = $2;
		$ranges{$lower}=$upper;
	}
}

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
	`java -Xmx2g -jar $gatk_path -R $ref -T DepthOfCoverage -o $read_dir$output -I $read_dir$file -omitIntervals -ct 1 -ct 10 -ct 30 -ct 50 -ct 75 -ct 100`;

	# Read results into Array and get rid of the first (header) line
	open OUT, '<', $read_dir.$output;
	my @lines = <OUT>;
	close OUT;
	shift @lines;

	my %base_of_depth;

	foreach my $line (@lines)
	{
		my @split = split(/\t/,$line);
		$split[0] =~ m/:(\d+)$/;
		my $pos = $1;
		my $depth = $split[1];
		$depth_arr[--$pos] = $depth;
		if($run_with_range) 
		{
			foreach my $lower (keys(%ranges))
			{
				my $upper = $ranges{$lower};
				next if($lower <= $pos && $pos >= $upper )
			}
		}
		$base_of_depth{$depth}++;
	}
	open POS, ">", $output_dir.$output.".position";

	my $counter = 1;	
	foreach my $dep (@depth_arr)
 	{
		say POS "$counter $dep";
 		$counter++;
	}
	close POS;
	
	open COV, ">", $output_dir.$output.".coverage";
	say COV "Depth\t$file\tCumulative";
	my $cumulative = 0;	
	foreach my $depth (sort { $a <=> $b } keys(%base_of_depth))
	{		
		$cumulative += $base_of_depth{$depth};
		say COV "$depth\t$base_of_depth{$depth}\t$cumulative";
	}
	close COV;
}

#Clean UP!
`rm $read_dir*Depth.*`;
