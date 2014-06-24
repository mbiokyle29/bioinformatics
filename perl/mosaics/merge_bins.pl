#!/usr/bin/env perl
# Merge bins of genome files into one huge one to feed to mosaics
# Takes a dir and a type (map, N, GC,)
# Expects 'mosaics' file like file names
use warnings;
use strict;
use Data::Printer;
use feature qw|say|;
use Getopt::Long;
use File::Slurp;

# Predef
my ($dir, $type, $output);
GetOptions(
	"dir=s" => \$dir,
	"type=s" => \$type,
	"output=s" => \$output,
);

unless($type =~ m/GC|N|map/ && $dir)
{
	die "merge_bins.pl --dir DIR --type (GC|N|map) --output output";
}

my @files = grep {/$type/} read_dir($dir);
p(@files);

open my $outfile, ">", $dir.$output;

foreach my $file (@files)
{
	$file =~ m/(^chr.+)(?=_$type)/;
	my $chrID = $1;
	say "$file generated id = $chrID";
	my @data = read_file($file);
	foreach my $dataline (@data)
	{
		chomp($dataline);
		say $outfile "$chrID\t$dataline";
	}
}
close $outfile;
