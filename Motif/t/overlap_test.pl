#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;
use lib "lib/";
use Peak;


my %cobound = (
	ebna2   => 0,
	ebna3a  => 0,
	ebna3b  => 0,
	ebna3c  => 0,
	rbpj92  => 0,
	rbpj234 => 0
);

my %cobound1 = (
	ebna2   => 0,
	ebna3a  => 0,
	ebna3b  => 0,
	ebna3c  => 0,
	rbpj92  => 0,
	rbpj234 => 0
);

my $peak_a = Peak->new(
	id         => "1",
	chromosome => "a",
	start      => 1,
	stop       => 10,
	summit     => 5,
	sequence   => "TTTTCCCCC",
	cobound    => \%cobound
);


my $peak_b = Peak->new(
	id         => "2",
	chromosome => "a",
	start      => 20,
	stop       => 30,
	summit     => 40,
	sequence   => "TTTAAAAAAA",
	cobound    => \%cobound1
);

my $res = $peak_a->overlaps($peak_b);
say $res;