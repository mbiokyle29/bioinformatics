#!/usr/bin/perl
use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Cigar;
use Data::Dumper;

my $string = shift;
my $cigar = Cigar->new(
	raw_string => $string,
	start_pos => 10
);

print Dumper @{ $cigar->array } ;
print $cigar->stack."\n";
print $cigar->length."\n";