#!/usr/bin/perl
# Check a merged bin file to make sure inc by 50 each
use warnings;
use strict;
use File::Slurp;

# Pass file and bin size

my $file = shift;
my $bin_size = shift;
my @lines = read_file($file);

my $last;
my @zero = split(" ", pop(@lines));

if($zero[1] == 0)
{
	$last = $zero[1];
} else { die "$file didnt start with zero!\n"; }

foreach my $line_scale (@lines)
{
	print $line_scale ."\n";
	my @line = split(" ", $line_scale);
	my $curr = $line[1];
	print "$last  $line[1]\n";
	if(($line[1] - $bin_size) == $last)
	{
		$last = $line[1];
	} else {
		die "$line[1] is misaligned!!\n";
	}
} 