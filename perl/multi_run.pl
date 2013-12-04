#!/usr/bin/perl
#
# Kyle McChesney
# Start of pipeline, run a folder full of reads through bowtie2
#
use warnings;
use strict;
use Getopt::Long;
use threads;
use Cwd qw(abs_path);

my $read_dir;
my $map_base;
my $matched = 0;

GetOptions ("read_dir=s" => \$read_dir,
			"map_base=s" => \$map_base,
			 "matched=i" => \$matched) or die("malformed command line args \n");


# Grab everything after bin name for fullpath
my $fp_bowtie2X = `which bowtie2`;
chomp($fp_bowtie2);
my $fp_results = $fp_bowtie2."results/";
my $fp_read = $fp_bowtie2."reads/";

print $fp_read."\n";
print $fp_bowtie2."\n";

opendir DIR, $fp_read;
my @read_files = grep { !/^\./ } readdir(DIR);
closedir DIR;

my $match_arg = "-U" unless $matched;

for my $read_file (@read_files)
{	
	if($read_file =~ m/^(.*)\.fastq$/)
	{
		my $output = $1.".sam";
		my $results = "bowtie2 -t --no-unal -x $map_base $fp_read$read_file -S $fp_results$output";
		print "For $read_file: \n";
		print $results."\n";
	}
}

