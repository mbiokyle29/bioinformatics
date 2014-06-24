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
# Dependencies
# -bowtie2
# -nproc or sysctil
#
# TODO
# Need to fix map full path, and output files, better implement skip if EBV genome
# Error messages and help also...

use warnings;
use strict;
use Getopt::Long;
use feature qw(say);
use Switch;

my $read_dir;
my $map_base;
my $matched = 0;
my $stat;

GetOptions(
    "d=s"       => \$read_dir,
    "m=s"       => \$map_base,
    "matched=i" => \$matched,
    "stats=i" => \$stat
) or die("malformed command line args \n");

# Figure out the number of cores to run on (total on machine - 1)
my $cores;
my $kernel = `uname -s`;
chomp($kernel);

switch ($kernel) {
    case "Linux"  { $cores = `nproc`; }
    case "Darwin" { $cores = `sysctl -n hw.ncpu`; }
    else          { die "Get a different Kernel \n"; }
}

chomp($cores);

# $genoms{'genomeX'} = (count of that genome in alignment)
our %genomes;

# Grab path for bowtie
my $fp_bowtie2 = "/home/kyle/lab/bowtie2/";
my $fp_results = $fp_bowtie2 . "results/";
my $fp_maps    = $fp_bowtie2 . "maps/";

# Grab all the files, ignore . and ..
my $dir_h;
opendir $dir_h, $read_dir;
my @read_files = grep { !/^\./ } readdir($dir_h);
closedir $dir_h;

# Default to unmatched unless specified in command
my $match_arg = "-U" unless $matched;

for my $read_file (@read_files) {
    my $read_count = 0;
    if ( $read_file =~ m/^(.*)\.f(ast)?q$/ ) {
        my $base_name = $1;
        my $fp_sam    = $fp_results . $1 . ".sam";
        my $results =
"bowtie2 -p $cores -t --no-unal -x $fp_maps$map_base $read_dir$read_file -S $fp_sam";
        my $bowtie_output = `$results`;
        say $bowtie_output;

        # Run Sam->Bam conversion
        &sam_to_bam( $fp_results . $base_name );

        # If aligning to EBV no need to run chromosome stats
        next unless($stat);

        # Output file and alignment counter
        my $ts          = time();
        my $sam_stat    = $fp_sam . ".$ts.stat";
        my $align_count = 0;

        # Read from one, write stats into other
        open my $sam_fh,  '<', $fp_sam;
        open my $stat_fh, '>', $sam_stat;

        # read in one line at a time and process
        while (<$sam_fh>) {
            my $in = $_;
            next if ( $in =~ /^@/ );
            chomp($in);

            my @line = split( "\t", $in );
            $align_count++;

            # increment or add to hash
            if ( $genomes{ $line[2] } ) {
                $genomes{ $line[2] }++;
            }
            else { $genomes{ $line[2] } = 1; }
        }
        close $sam_fh;

        # Write the .stat file
        say $stat_fh "Chromosome frequencey results for $fp_sam:";
        for my $genome ( keys(%genomes) ) {
            say $stat_fh "$genome had $genomes{$genome} reads";
        }
        say $stat_fh "total alignements = $align_count";
        say $stat_fh "output from bowtie2: ";
        say $stat_fh "$bowtie_output";
        close $stat_fh;
    }
}

# Takes filename and full path, makes bam files
# TODO should prolly handle errors and return codes
sub sam_to_bam {
    my $base_name = shift;
    my $sam       = $base_name . ".sam";
    my $bam       = $base_name . ".bam";
    my $sorted    = $base_name . ".sorted";
    `samtools view -S -b -o $bam $sam`;

    # That speed up though
    `samtools-rs rocksort -@ 8 -m 16G $bam $sorted`;
    `samtools index $sorted.bam`;
}
