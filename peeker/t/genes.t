#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use Test::Simple tests => 2;
use Test::Exception;
use lib "../lib";
use Gene;

my $gene = Gene->new(gene_start => 1, gene_stop => 2, name => "GTD", gtf => "hg19");
ok(defined($gene) && ref($gene) eq "Gene", "Gene init okay");

$gene = Gene->new(gene_start => 10, gene_stop => 2, name => "GTD", gtf => "hg19");
ok(defined($gene) && ref($gene) eq "Gene", "Reverse strand gene created okay");

dies_ok {
	my $gene = Gene->new(gene_start => 1, gene_stop => -2, name => "GTD", gtf => "hg19");	
} "Died with negative stop";
