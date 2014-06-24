#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use feature qw|say|;

# Example: bowsaics.pl -d /dir/with/read/ -m EBV -i input -c ebna2chip
my ($read_dir, $map_header, $chip, $input);
GetOptions (
	# Directory with the raw reads (fastq fa fq etc)
	"d=s" => \$read_dir, 

	# Name of the bowtie2 map
	"m=s" => \$map_header,
	
	# Basename (not .fa) of input file
	"i=s" => \$input,
	
	# Basename (not .fa) of chip file
	"c=s" => \$chip
);

my $fp_bowtie2 = "/home/kyle/lab/bowtie2/";
my $fp_results = $fp_bowtie2 . "results/";
$chip .= ".sam";
$input .= ".sam";

# Run the bowtie script
say `perl /home/kyle/lab/perlpipe/perl/bowtie/multi_run.pl -d $read_dir -m $map_header`;

# Move the sams from the bowtie folder to here
move($fp_results.$chip, "./");
move($fp_results.$input, "./");

# Run mosaics
say `perl /home/kyle/lab/perlpipe/perl/mosaics/MosaicsPipe.pl --type IO --format sam --chip $chip --input $input`;
