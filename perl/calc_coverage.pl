#!/usr/bin/perl
# This script takes a full file path to a directory with alignment results in
# sam format. It sorts through each file, and builds a hash as follows:
# $base_counts{ # A given position on the chromosome } = # The number of reads that contain a base at this position
# It then outputs the results in an organized fashion and creates a wig file of the results
#
# NOTE! SAM results must all be of the same chromosome!
# Kyle McChesney

use warnings;
use strict;
use feature qw(say);
use Data::Dumper;

# Directory with sam files
my $read_dir = shift;

# Array to build
## NOTE THIS IS HARD CODED
my $total = 171823;
my %base_counts;

for(my $i = 1; $i <= $total; $i++)
{
	$base_counts{$i} = 0;
}

# Name of chromosome in SAMs, found during processing
my $chromosome;

# Get only the .sam files
opendir DIR, $read_dir;
my @align_files = grep { /.bam$/ } readdir(DIR);
closedir DIR;

foreach my $file (@align_files)
{
	my $
}