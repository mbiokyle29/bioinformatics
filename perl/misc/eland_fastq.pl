#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use feature qw(say);

my $eland_file = shift;
#my @lines = read_file($eland_file);
open my $file, "<", $eland_file;
while(<$file>)
{
	my $line = $_;
	chomp($line);
	my @line_arr = split(/\s+/, $line);
	if(scalar(@line_arr) == 4)
	{
		next if ($line_arr[2] =~ m/QC/);
		say "@".$line_arr[0];
		say $line_arr[1];
		say "+".$line_arr[0];
		say "NA";
	}
	else
	{
		next if ($line_arr[10] =~ m/QC/);
		say "@".$line_arr[0];
		say $line_arr[8];
		say "+".$line_arr[0];
		say $line_arr[9];
	}
}