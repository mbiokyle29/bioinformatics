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

# $genoms{'genomeX'} = (count of that genome in alignment)
our %genomes;

GetOptions ("read_dir=s" => \$read_dir,
			"map_base=s" => \$map_base,
			 "matched=i" => \$matched) or die("malformed command line args \n");

# Grab path for bowtie
my $fp_bowtie2 = `which bowtie2`;
chomp($fp_bowtie2);

# Build full paths 
# needs reads and results folder in bowtie2 dir
my $fp_results = (substr $fp_bowtie2, 0, -length("bowtie2"))."results/";
my $fp_read = (substr $fp_bowtie2, 0, -length("bowtie2"))."reads/".$read_dir."/";

# Figure out the number of cores to run on (total on machine - 1)
my $cores = `nproc` - 1; 

# Grab all the files, ignore . and ..
opendir DIR, $fp_read;
my @read_files = grep { !/^\./ } readdir(DIR);
closedir DIR;

# Default to unmatched unless specified in command
my $match_arg = "-U" unless $matched;

for my $read_file (@read_files)
{	
	my $read_count = 0;
	if($read_file =~ m/^(.*)\.fastq$/)
	{
		my $output = $1.".sam";
		my $results = "bowtie2 -p $cores -t --no-unal -x $map_base $fp_read$read_file -S $fp_results$output";
		print "For $read_file: \n";
		print $results."\n";
		
		# Output file and alignment counter
		my $sam_stat = $output.".stat";
		my $align_count = 0;
		
		# Read from one, write stats into other
		open  SAM, '<', $output;
		open  STAT, '>', $sam_stat;
		
		# one line at a time
		while (<SAM>)
		{
			my @line = split("\t", chomp($_));
			if exists ($genomes{$line[2]})
			{
				$genomes{$line[2]}++;
			}
			else
			{
				$genomes{$line[2]} = 1;
			}
		}
	}
}

























