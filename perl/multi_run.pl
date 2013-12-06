#!/usr/bin/perl
#
# Kyle McChesney
# Start of pipeline, run a folder full of reads through bowtie2
# need bowtie2 and nproc installed
use warnings;
use strict;
use Getopt::Long;
use threads;
use feature 'say';

my $read_dir;
my $map_base;
my $matched = 0;

# $genoms{'genomeX'} = (count of that genome in alignment)
our %genomes;

GetOptions ("d=s" => \$read_dir,
			"m=s" => \$map_base,
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
		
		# read in one line at a time and process
		while (<SAM>)
		{
			my @line = split("\t", chomp($_));
			$align_count++;
			# increment or add to hash
			if ($genomes{$line[2]})
			{
				$genomes{$line[2]}++;
			} else { $genomes{$line[2]} = 1; }
		}
		close SAM;
		
		# Write the .stat file
		say STAT "Chromosome frequencey results for $output:";
		for my $genome (keys(%genomes))
		{
			say STAT "$genome had $genomes{$genome} reads";
		}
		say STAT "total alignements = $align_count";
		close STAT;
		
		
	}
}

























