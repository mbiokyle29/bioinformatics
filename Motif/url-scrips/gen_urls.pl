#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;
use List::Util qw|min max|;
use lib "/home/kyle/lab/Binder/lib";
use Binder;

my @misfits = read_file(shift);
my $base= "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=(?)%3A(?)-(?)&hgsid=384337607_fNiMyItEp1Bz0ZVRt1A9NBASaHzY";
my $binder = Binder->new( base_string => $base);
shift(@misfits);

my @starts;
my @stops;
my $chr;

while(my $line = shift(@misfits))
{
	if($line eq "\n")
	{
		shift(@misfits);
		my $min = min(@starts);
		my $max = max(@stops);

		$min -= 100;
		$max += 100;

		$binder->bind($chr, $min, $max);
		say $binder->bound_string;
		undef(@starts);
		undef(@stops);
	}

	else
	{
		my ($t_chr, $t_start, $t_stop) = (split("\t", $line))[0..2];
		push(@starts, $t_start);
		push(@stops, $t_stop);
		$chr = $t_chr;
	}
}