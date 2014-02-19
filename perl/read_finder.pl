#!/usr/bin/perl
#90985 90986 90987
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;

my ($start_pos, $end_pos, $sam);

GetOptions (
	"s=i" => \$start_pos,
	"e=i" => \$end_pos,
	"r=s" => \$sam,
);

# CHECK ARGS
unless ($start_pos && $end_pos && $sam)
{
	die "read_finder.pl -s starting_position -e ending_position -r sam_files";
}

if($end_pos < $start_pos)
{
	die "ending position cannot be less then start position!";
}

# Get reads
my @lines = read_file($sam);
foreach my $line (@lines)
{	
	next if ($line =~ m/^@/);
	my @fields = split("\t", $line);

	my $seq = $fields[9];
	my $read_start = $fields[3];
	my $read_length = (length($seq)-1);
	my $read_end =  $read_start+$read_length;
	
	# now smart match
}


