#!/usr/bin/perl
use warnings;
use strict;
use Statistics::R;

my $sam = shift;
my $r_con = Statistics::R->new();
$r_con->run("library(mosaics)");
my $wiggle_command = "generateWig( infile=\"$sam\", fileFormat=\"sam\", outfileLoc=\"./\")";
$r_con->run($wiggle_command);