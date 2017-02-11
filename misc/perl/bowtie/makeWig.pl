#!/usr/bin/env perl
#
# Kyle McChesney
# testTform.t
# 
#
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use feature 'say';
use FileDo;

my $sam = shift;
$sam =~ m/(.*)\.sam/;
my $basename = $1;

my $size_file = shift;
$size_file =~ m/.*\/(.*)\.size/;
my $chrom = $1;

my $sorted_bam = sam_to_bam($sam,$basename);
my $bed = $basename.".bed";
my $wig = $basename.".wig";

`bamToBed -i $sorted_bam > $bed`;
`genomeCoverageBed -d -i $bed -g $size_file >> $wig`;
`awk '{print \$3}' $wig > $wig.temp`;
my $data = read_file($wig.".temp");

my $track_line = "track type=wiggle_0 name=\"(?) \" description=\"(?) wig file\"\n";
$track_line .= "fixedStep chrom=$1 start=1 step=1\n";
my $binder = Binder->new(base_string => $track_line);
$binder->bind($sam,$sam);
write_file($wig.".clean", $binder->bound_string());
append_file($wig.".clean", $data);
`rm $wig $wig.temp $bed`;

# Takes filename and full path, makes bam files
# TODO should prolly handle errors and return codes
sub sam_to_bam {
    my $sam       = shift;
   	my $basename = shift;
    my $bam       = $basename.".bam";
    my $sorted    = $basename.".sorted";
    `samtools view -S -b -o $bam $sam`;

    # That speed up though
    `samtools-rs rocksort -@ 8 -m 16G $bam $sorted`;
    say "Returning $sorted as the sorted bam file";
    return $sorted.".bam";
}
