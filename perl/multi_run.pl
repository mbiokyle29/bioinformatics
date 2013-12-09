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
use feature qw(say);
use Switch;

my $read_dir;
my $map_base;
my $matched = 0;

# Figure out the number of cores to run on (total on machine - 1)
my $cores;
my $os = `uname -s`;
chomp($os);

switch($os)
{
	case "Linux" { $cores = `nproc`; }
	case "Darwin" { $cores = `sysctl -n hw.ncpu`; }
	else die "Are you on windows? \n";
}

chomp($cores);

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
my $fp = (substr $fp_bowtie2, 0, -length("bowtie2"));
my $fp_results = $fp."results/";
my $fp_read = $fp."reads/".$read_dir."/";
my $fp_maps = $fp."maps/";



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
		my $base_name = $1;		
		my $fp_sam = $fp_results.$1.".sam";
		my $results = "bowtie2 -p $cores -t --no-unal -x $fp_maps$map_base $fp_read$read_file -S $fp_results$sam";
		my $bowtie_output = `$results`;
		say $bowtie_output;
		
		# Run Sam->Bam conversion
		&sam_to_bam($fp_results.$base_name);

		# If aligning to EBV no need to run chromosome stats
		next if($map_base eq "EBV");
		
		# Output file and alignment counter
		my $ts = time();
		my $sam_stat = $fp_$sam.".$ts.stat";
		my $align_count = 0;
		
		# Read from one, write stats into other
		open  SAM, '<', $fp_sam;
		open  STAT, '>', $fp_stat;
		
		# read in one line at a time and process
		while (<SAM>)
		{
			my $in = $_;
			next if ($in =~ /^@/);
			chomp($in);
			
			my @line = split("\t", $in);
			$align_count++;
			# increment or add to hash
			if ($genomes{$line[2]})
			{
				$genomes{$line[2]}++;
			} else { $genomes{$line[2]} = 1; }
		}
		close SAM;
		
		# Write the .stat file
		say STAT "Chromosome frequencey results for $sam:";
		for my $genome (keys(%genomes))
		{
			say STAT "$genome had $genomes{$genome} reads";
		}
		say STAT "total alignements = $align_count";
		say STAT "output from bowtie2: ";
		say STAT "$bowtie_output";
		close STAT;
	}
	
	# Takes filename and full path, makes bam files
	# TODO should prolly handle errors and return codes
	sub sam_to_bam 
	{
		my $base_name = shift;
		my $sam = $base_name.".sam";
		my $bam = $base_name.".bam";
		my $sorted = $base_name.".sorted";
		`samtools view -S -b -o $bam $sam`;
		`samtools sort $bam $sorted`;
		`samtools index $bam`;
	}
}

