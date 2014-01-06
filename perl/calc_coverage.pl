#!/usr/bin/perl
# This script takes a full file path to a directory with alignment results in
# sam format. It sorts through each file, and builds a hash as follows:
# $base_counts{ # A given position on the chromosome } = # The number of reads that contain a base at this position
# It then outputs the results in an organized fashion
# Kyle McChesney

use warnings;
use strict;
use Getopt::Long;
use feature qw(say);
use Data::Dumper;

my $read_dir = shift;

my %base_counts;

opendir DIR, $read_dir;
my @align_files = grep { !/^\./ } readdir(DIR);
closedir DIR;

foreach my $file (@align_files)
{
	next if($file !~ m/\.sam$/);
	open SAM, '<', $read_dir.$file;
	while(<SAM>)
	{
		next if $_ =~ m/^@/;
		my @line = split(/\t/,$_);
		my $start_pos = $line[3];
		my $read = $line[9];
		chomp($read);
		my $length = length($read) - 1;
		my $counter = $start_pos + $length;
		while($counter >= $start_pos)
		{
			$base_counts{$counter}++;
			$counter--;
		}
	}
}

say "Poisition   # of reads";
foreach my $key (sort {$a<=>$b} keys(%base_counts))
{
	say "$key            $base_counts{$key}";
}