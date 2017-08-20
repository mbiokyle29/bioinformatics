#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;

my @urls = read_file(shift);
my $out = "peaks_merged.html";
my $out_stuff;

$out_stuff .= "<!DOCTYPE html>\n";
$out_stuff .= "<html lang=\"en-US\">\n";
$out_stuff .= "<head></head>";
$out_stuff .= "<body>\n\t<h1>Peak Links</h1>\n";
my $counter = 1;
foreach my $url (@urls)
{
	chomp($url);
	$out_stuff .= "\t<a href=\"$url\">peak #$counter</a><br>\n";
	$counter++;
}
$out_stuff .= "</body>";

write_file($out, $out_stuff);