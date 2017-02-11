#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;

my $sam_file = shift;
my @lines = read_file($sam_file);
my $output = "out.fa";

foreach my $line (grep(!/^@/, @lines))
{
	my @sam_fields = split("\t", $line);
	my $sam_seq = $sam_fields[9];
	my $outstring = ">fa conversion from $sam_file";
	$outstring .= $sam_seq."\n";
	append_file($output, $outstring);
}