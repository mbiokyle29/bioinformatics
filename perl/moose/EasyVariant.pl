#!/usr/bin/perl
use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Genome;
use Cigar;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say switch);

my ($start_pos, $end_pos, $sam, $reference);

GetOptions (
	"start=i" => \$start_pos,
	"end=i" => \$end_pos,
	"sam=s" => \$sam,
	"reference=s" => \$reference,
	"bad=i" => \$keep_bad
);

my $genome = Genome->new(
	fasta => $reference
);

my $cigar = Cigar->new(
	raw_string => "3M",
);

say $genome->length;
say $genome->seq;
say $genome->base_at(1);
say $genome->base_at(2);
say $genome->base_at(3);
say $genome->base_at(4);
say $genome->base_at(5);
say $genome->base_at(6);
say $genome->base_at(10);
say $genome->base_at(10);
say $genome->base_at(10);
say $genome->base_at(10);