#!/usr/bin/perl
use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Cigar;
use Genome;
use Read;
use Data::Dumper;

my $genome = Genome->new(
	fasta => shift,
);

my $cigar = Cigar->new(
	raw_string => "3M",
);

print $genome->length;