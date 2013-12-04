#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use threads;

my $read_dir;
my $map_base;
my $matched = 0;

GetOptions ("read_dir=s" => \$read_dir,
			"map_base=s" => \$map_base,
			 "matched=i" => \$matched) or die("malformed command line args \n");

opendir DIR, $read_dir;
my @read_files = grep { !/^\./ } readdir(DIR);
closedir DIR;

my $match_arg = "-U" unless $matched;

for my $read_file (@read_files)
{	

	if($read_file =~ m/^(.*)\.fastq$/)
	{
		my $output = $1.".sam";
		my $results = "bowtie2 -t --no-unal -x $map_base $read_dir$read_file -S $output";
		print "For $read_file: \n";
		print $results."\n";
	}
}

