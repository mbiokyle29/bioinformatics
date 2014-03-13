#!/usr/bin/perl
use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Cigar;

my $string = shift;
my $cigar = Cigar->new(
	raw_string => $string
);
print $cigar->stack;
print $cigar->length;