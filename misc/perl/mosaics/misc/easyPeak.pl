#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use feature qw|say|;
use Data::Printer;

my $chip = shift;
my $input = shift;
my $slide = shift;
my $pointer = 0;

my @chips = map($_[2], split("\t", read_file($chip)));

my @inputs = read_file($input);

my $chip_length = @chips;
my $input_length = @inputs;
my $shortest = $chip_length > $input_length ? $chip_length : $input_length;
my @chip_chunk;
my @input_chunk;

while(($pointer + $slide) > $chip_length)
{
	@chip_chunk = (@chips[$pointer .. $pointer+$slide]);
	@input_chunk = (@chips[$pointer .. $pointer+$slide]);
	$pointer += $slide;

	my ($chip_total, $input_total) = 0;
	

}