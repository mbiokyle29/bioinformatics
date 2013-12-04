#!/usr/bin/perl
use warnings;
use strict;

my $dir = shift;
my $result = shift;

open FILE, "> $result";
close FILE;

opendir DIR, $dir;
my @files = grep { !/^\.+/ } readdir DIR;

for my $file (@files)
{
	print "going to cat $file\n";
	`grep -v '>' $dir/$file >> $result`;
}
