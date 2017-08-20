#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;
use DBI;
use Getopt::Long;

# DBI connect
my $dbh = DBI->connect('DBI:mysql:Peaks', 'peaker', '11peaks12');
my @proteins;
my @nots;
my @ors;
my $trim;

GetOptions(
	"with=s" => \@proteins,
	"or=s" => \@ors,
	"not=s" => \@nots,
	"trim=i" => \$trim
);

# Base Query
my $query_base = "SELECT * FROM refined_peaks ";
my ($and_addition, $or_addition, $not_addition);

# Handle ANDs
if(@proteins)
{
	foreach my $prot (@proteins)
	{
		unless(&valid_prot($prot)) { die "$prot is not a valid protein!"; }

		unless($and_addition) 
		{ 
			$and_addition = "( $prot=1 "; 
		} else { 
			$and_addition .= "AND $prot=1 "; 
		}
	}
	$and_addition .= ") ";
}

# Handle (ORs)
if(@ors)
{
	foreach my $or (@ors) {

		unless(&valid_prot($or)) { die "$or is not a valid protein!"; }
		
		unless($or_addition) {
			$or_addition = "($or=1 ";
		} else {
			$or_addition .= "OR $or=1 ";
		}
	}
	$or_addition .= ") ";
}

# Handle not
if(@nots)
{
	foreach my $not (@nots)
	{
		unless(&valid_prot($not)) { die "$not is not a valid protein!"; }

		unless($not_addition)
		{
			$not_addition = "( $not=0 ";	
		} else {
			$not_addition .= "AND $not=0 ";
		}
	}
	$not_addition .= ") ";
}

if($and_addition or $or_addition or $not_addition) {
	$query_base .= "WHERE ";
	
	if($and_addition) {
		$query_base .= $and_addition." AND ";
	}

	if($or_addition) {
		$query_base .= $or_addition." AND ";
	}

	if($not_addition) {
		$query_base .= $not_addition;
	}
}

$query_base =~ s/\s+AND\s+$//g;
say $query_base;
die;

my $file = "peaks";
$file .= "_".$_ for (@proteins);
$file .= "_or_".$_ for (@ors);
$file .= "_not-".$_ for (@nots);
$file .= ".fa";

my $sth = $dbh->prepare($query_base);
$sth->execute();
while(my $row = $sth->fetchrow_hashref())
{
	my $fasta_header = ">";
	my $start = $$row{start};
	my $stop = $$row{stop};

	$fasta_header .= $$row{chromosome}."|";
	$fasta_header.= $start."-".$stop."|";

	my $cobound_id = "";
	$cobound_id .= $_."-" for( grep{$$row{$_} == 1} @proteins);
	chop($cobound_id);

	my $sequence = $$row{sequence};
	
	if($trim)
	{
		$sequence = trim_seq($sequence, $trim);
		$start += $trim;
		$stop -= $trim;
	}

	$fasta_header .= $cobound_id."|\n";
	$fasta_header .= $sequence."\n";
	append_file($file, $fasta_header);
}

my $meme_command = "meme-chip -o $file-OUT/ -meme-p 15 ";
$meme_command .= "-db /home/kyle/Applications/meme/db/JASPAR_CORE_2014_vertebrates.meme ";
$meme_command .= "-db /home/kyle/Applications/meme/db/uniprobe_mouse.meme ";
$meme_command .= "$file";
say `$meme_command`;

sub trim_seq {
	my ($seq, $trim) = @_;
	my $neg = $trim * (-1);
	return(substr($seq, $trim, $neg));
}

sub valid_prot {
	my ($prot) = @_;
	my @valid_proteins = qw|ebna3a ebna2 ebna3b ebna3c rbpj92 rbpj234|;
	return ($prot ~~ @valid_proteins);
}