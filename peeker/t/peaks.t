#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use Test::Simple tests => 5;
use Test::Exception;
use lib "../lib";
use Peak;

my $peak = Peak->new(
	peak_start => 1, 
	peak_stop => 3,
	gtf => "hg19",
	file => "peaks.pk"
);

ok(defined($peak) && ref($peak) eq "Peak", "Peak init okay");
ok($peak->peak_summit == 2, "Peak summit set to average correcty");

$peak = Peak->new(
	peak_start => 1, 
	peak_stop => 3,
	peak_summit => 2,
	gtf => "hg19",
	file => "peaks.pk"
);

ok($peak->peak_summit == 2, "Peak summit set by constructor correcty");

dies_ok {
	my $bad_peak = Peak->new(
		peak_start => 5, 
		peak_stop => 3,
		peak_summit => 2,
		gtf => "hg19",
		file => "peaks.pk"
	)
} "Failed correctly with bad peak start";

dies_ok {
	my $bad_peak = Peak->new(
		peak_start => 1, 
		peak_stop => 3,
		peak_summit => 10,
		gtf => "hg19",
		file => "peaks.pk"
	)
} "Failed correctly with bad peak summit";