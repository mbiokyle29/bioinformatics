#!/usr/bin/perl
#
# Kyle McChesney
# Start of pipeline, run a folder full of reads through bowtie2
# need bowtie2 and nproc installed
#
# Need to fix map full path, and output files, better implement skip if EBV genome
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
# for linux my $cores = `nproc` - 1; 
my $cores = `sysctl -n hw.ncpu`;
chomp($cores);

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
		my $bowtie_output = `$results`;
		say $bowtie_output;

		next if($map_base eq "EBV");
		# Output file and alignment counter
		my $ts = time();
		my $sam_stat = $fp_results.$output.".$ts.stat";
		my $align_count = 0;
		my $fp_output = $fp_results.$output;
	
		# Read from one, write stats into other
		open  SAM, '<', $fp_output;
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
		say STAT "output from bowtie2: ";
		say STAT "$bowtie_output";
		close STAT;
	}
}

