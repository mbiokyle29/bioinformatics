#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;
use Getopt::Long;
use Genome;
use Sequence;
use constant PROMOTER_LENGTH => 7;
my ($fasta, $before, $promoter_file);

# Fasta file, Promoter file, amount before promoter to grab
# gene file is tsv: name<tab>cds-start

GetOptions(
	"fasta=s"  => \$fasta,
	"promoters=s"  => \$promoter_file,
	"before=i" => \$before,
	#"help=?"   => \&useage(),
);

# Verify
unless ($fasta && $promoter_file && $before) {
	&useage("--fasta --promoters --before flags are required");
}

# get genome file into mysql
my $genome = Genome->new(fasta => $fasta);
say $genome->name();
my @promoters = read_file($promoter_file);

# split the list into forward and rev promoters
my @for_promoters = grep { /RF/i } @promoters;
my @rev_promoters = grep { /LF/i } @promoters;

# Do forward first
# want range (Start - Before) --> Start
# note we dont want the actual base of the promoter
foreach my $gene (@for_promoters) {
	my ($name, $prom) = split(/\t/, $gene);
	my $start = ($prom - $before);
	my $seq_raw = get_geneomic_region($start, ($prom-1), $genome);
	my $seq = Sequence->new(seq => $seq_raw, name => $name);
	my $fasta = ">".$seq->name()."- $before bases before the promoter\n";
	$fasta.= $seq->seq();
	write_file($seq->name().".fasta", $fasta);
}

# Do then do rev
# same idea except grab the sequence after
# and get the rev complement
foreach my $gene (@rev_promoters) {
	my ($name, $prom) = split(/\t/, $gene);
	my $end = ($prom+PROMOTER_LENGTH + $before);
	my $seq_raw = get_geneomic_region(($prom+PROMOTER_LENGTH+1), $end, $genome);
	$seq_raw = rev_complement($seq_raw);
	my $seq = Sequence->new(seq => $seq_raw, name => $name);
	my $fasta = ">".$seq->name()."- $before bases before the promoter REV COM\n";
	$fasta.= $seq->seq();
	write_file($seq->name().".fasta", $fasta);
}

# subroutine for generating the sequences
sub get_geneomic_region {
	my ($start, $stop, $genome) = @_;
	$genome->sequence_range($start,$stop);
}

# Compute reverse complement of DNA strand
# using tr///
sub rev_complement {
	my $seq = shift;
	my $rev = reverse($seq);
	$rev =~ tr/ACGTacgt/TGCAtgca/;
	return $rev;
}