#!/usr/bin/perl
#
# Kyle McChesney
# Start of a bowtie alignment pipeline
# Also does so simple chromosome analysis when aligning against multi-chromosome maps
# Is picky about file structure
# -takes a directory (fullpath) full of read files in fastq format
# -takes a map_base which is the base name of the bowtie mapped reference to use
# -can also set bowtie2 to run in matched mode with -matched flag
#
# TODO
# Need to fix map full path, and output files, better implement skip if EBV genome
# Error messages and help also...
use warnings;
use strict;
use Getopt::Long;
use feature qw(say);
use Switch;
use File::Slurp;

# Directories and references
my ($read_dir, $map_dir, $map_name, $fasta_file, $results_dir, $script_path);

# Options
my ($cores, $matched);

# Optional load from file
my $config_file;

GetOptions(
    "read-dir=s"    => \$read_dir,
    "map-dir=s"     => \$map_dir,
    "map-name=s"    => \$map_name,
    "results-dir=s" => \$results_dir,
    "cores=i"       => \$cores,
    "matched=i"     => \$matched,
    "script-path=s"   => \$script_path,
    "config=s"      => \$config_file,
    "fasta=s"      => \$fasta_file
) or die("malformed command line args \n");

# Do config file
if($config_file) {
    validate_config_file($config_file);
}

# Grab all the files, ignore . and ..
my @read_files = read_dir($read_dir);

# Default to unmatched unless specified in command
my $match_arg = "-U" unless $matched;

for my $read_file (@read_files) {
    my $read_count = 0;
    if ( $read_file =~ m/^(.*)\.f(ast)?q(\.gz)?$/ ) {
        my $base_name = $1;
        my $fp_sam    = $results_dir . $1 . ".sam";
        my $results ="bowtie2 -p $cores -t --no-unal -x $map_dir$map_name $read_dir$read_file -S $fp_sam";
        my $bowtie_output = `$results`;
        say $bowtie_output;
    } else {
        say "Skipping non read file: $read_file";
        next;
    }
}

# Takes filename and full path, makes bam files
# TODO should prolly handle errors and return codes
sub sam_to_bam {
    my $base_name = shift;
    my $sam       = $base_name . ".sam";
    my $bam       = $base_name . ".bam";
    my $sorted    = $base_name . ".sorted.bam";
    `samtools view -S -b -o $bam $sam`;

    # That speed up though
    `samtools-rs rocksort -@ 8 -m 16G $bam $sorted`;
    `samtools index $sorted`;
    say "Returning $sorted as the sorted bam file";
    return $sorted;
}


sub useage {
    say "Bowtie2Pipe.pl useage: ";
    say "Bowtie2Pipe.pl -d DIRECTORY_WITH_READS -m BOWTIE_2_MAP_NAME --matched MATCHED_DATA_FLAG --stats STATISTICS_DATA_FLAG";
    say "-d and -m are required!"
}