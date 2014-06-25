#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use feature qw|say|;
use Data::Printer;

my $wiggle_file = shift;
my @wiggle;

if(-e $wiggle_file) {
 	@wiggle = read_file($wiggle_file) or die "Could not read $wiggle_file";
} else { die "File $wiggle_file does not exist"; }

my $track_header = shift(@wiggle);
my $buffer;
my $first = 1;
for my $line (@wiggle)
{
	if($line =~ m/^[a-z]+/i)
	{
		if($first)
		{
			$first = 0;	
		} else {
			&flush_buffer($buffer, $track_header);
		}
		$buffer = $line;
	} else {
		$buffer .= $line;
	}
}

&flush_buffer($buffer, $track_header);

sub flush_buffer {
	my ($buffer, $track_header) = @_;
	my $chr_line = ($buffer =~ /\A(.*?)$/ms)[0];
	$chr_line =~ m/\schrom=(\w+)\s/;
	my $chr_id = $1;
	$track_header =~ s/name=\"([^\"]+)\"/name=\"$1 split for $chr_id\"/;
	$track_header =~ s/description=\"([^\"]+)\"/description=\"$1 split for $chr_id\"/;
	$buffer = $track_header.$buffer;
	write_file($chr_id."splitWig.wig", $buffer);
}