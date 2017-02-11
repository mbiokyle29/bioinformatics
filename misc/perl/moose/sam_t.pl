#!/usr/bin/perl
use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Sam;
use Data::Dumper;
use File::Slurp;

my $sam = shift;
my $line = read_file($sam);
my $test = Sam->new(
	raw_string => $line
);

print Dumper $test;