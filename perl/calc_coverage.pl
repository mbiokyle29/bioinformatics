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

# Directory with sam files
my $read_dir = shift;
# Hash to build
my %base_counts;
# Name of chromosome in SAMs, found during processing
my $chromosome;

# Get only the .sam files
opendir DIR, $read_dir;
my @align_files = grep { /.sam$/ } readdir(DIR);
closedir DIR;

foreach my $file (@align_files)
{
	open SAM, '<', $read_dir.$file;
	while(<SAM>)
	{
		next if $_ =~ m/^@/;
		my @line = split(/\t/,$_);
		
		# Get the name of the chromosome once!
		unless($chromosome)
		{
			$chromosome = $line[2];
		}
		
		# Get the necessary fields
		my $start_pos = $line[3];
		my $read = $line[9];
		chomp($read);
		
		# Add length to start pos, count back down to start pos, incrementing
		my $length = length($read) - 1;
		my $counter = $start_pos + $length;
		while($counter >= $start_pos)
		{
			$base_counts{$counter}++;
			$counter--;
		}
	}
}

# Set up WIG file
my $ts = time();
my $wigfile = $read_dir.$chromosome.$ts."-wiggle.wig";
open WIG, '>', $wigfile;
say WIG "variableStep chrom=$chromosome";
say "Chromosome: $chromosome";
say "Poisition   # of reads";

# Hard coded percent coverage
my $total = 171823;
my $mapped = keys(%base_counts);
my $percent = ($mapped/$total)*100;

# Do both wig and output in same loop
foreach my $key (sort {$a<=>$b} keys(%base_counts))
{
	# Print to stdout
	say "$key            $base_counts{$key}";
	# Build wig file
	say WIG "$key $base_counts{$key}";
}
close WIG;
say "Percent mapped = % $percent";