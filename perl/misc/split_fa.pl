#!/usr/bin/perl
# Split a large read file into smaller chunk files
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;

my $in_file = shift;
my @big_file = read_file($in_file);
$in_file =~ s/\.[^.]+$//;

my $read_in = 0;
my $splits = 0;
my $buffer;

foreach my $line (@big_file)
{
	say "line";
	if($read_in != 200000)
	{
		$buffer .= $line;
		$read_in ++;
	}
	else {
		$splits ++;
		$read_in = 0;
		write_file($in_file."split-$splits", $buffer);
		$buffer = "";
	}
}
write_file($in_file."split-$splits", $buffer) if $buffer;