#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use feature qw|say|;
use Data::Printer;

my $bed_file = shift;
my @bed;

if(-e $bed_file) {
 	@bed = read_file($bed_file) or die "Could not read $bed_file";
} else { die "File $bed_file does not exist"; }

my $track_header = shift(@bed);
my ($buffer, $current_id, $line_id);
for my $line (@bed)
{
	$line_id = (split(/\s/, $line))[0];
	if(not $current_id) { 
		$current_id = $line_id;
		$buffer = $line;
	} elsif($line_id ne $current_id) {
		&flush_buffer($buffer, $track_header, $current_id);
		$current_id = $line_id;
		$buffer = $line;
	} else {
		$buffer .= $line;
	}
}

&flush_buffer($buffer, $track_header, $current_id);

sub flush_buffer {
	my ($buffer, $track_header, $chr_id) = @_;
	$track_header =~ s/name=([\S]+)/name=\"$1 split for $chr_id\"/;
	$track_header =~ s/description=\"([^\"]+)\"/description=\"$1 split for $chr_id\"/;
	$buffer = $track_header.$buffer;
	write_file($chr_id."splitBed.bed", $buffer);
}