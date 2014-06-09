#!/usr/bin/env perl
use warnings;
use strict;
use File::Slurp;
use Getopt::Long;
use Data::Printer;

unless(scalar(@ARGV) == 2)
{
	die "Please enter two files to compare";
}

# Read Files and dump header
my @file_one = read_file($ARGV[0]) or die "Could not read $ARGV[1]";
my @file_two = read_file($ARGV[1]) or die "Could not read $ARGV[2]";
shift(@file_one);
shift(@file_two);

# Merge
my @peaks = (@file_one, @file_two);

my %start_positions;
my %end_positions;

foreach my $line (@peaks)
{
	my @fields = split("\t", $line);
	$start_positions{$fields[0]}{$fields[1]}++;
	$end_positions{$fields[0]}{$fields[2]}++;
	
	if($start_positions{$fields[0]}{$fields[1]} == 2)
	{
		print "Matching peak start pos found: \n";
		print $fields[0]."\t".$fields[1]."\n";
	}

	if($end_positions{$fields[0]}{$fields[2]} == 2)
	{
		print "Matching peak end pos found: \n";
		print $fields[0]."\t".$fields[2]."\n";
	}
}
