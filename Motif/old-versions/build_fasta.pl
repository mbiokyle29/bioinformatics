#!/usr/bin/perl
use warnings; use strict;
use feature qw|say|;
use Data::Printer;
use DBI;
use File::Slurp;

# DBI connect
my $dbh = DBI->connect('DBI:mysql:Peaks', 'peaker', '11peaks12');

# Get all the refined peaks
my $get = "SELECT * FROM refined_peaks";
my $sth = $dbh->prepare($get);
$sth->execute();

my @proteins = qw|ebna3a ebna2 ebna3b ebna3c rbpj92 rbpj234|;
my @trims = qw|50 100 150|;
while(my $row = $sth->fetchrow_hashref())
{
	my $fasta_header = ">";
	my $start = $$row{start};
	my $stop = $$row{stop};
	$fasta_header .= $$row{chromosome}."|";
	my $no_trim = $fasta_header.$start."-".$stop."|";

	my $cobound_id = "";
	$cobound_id .= $_."-" for( grep{$$row{$_} == 1} @proteins);
	chop($cobound_id);

	my $sequence = $$row{sequence};
	$no_trim .= $cobound_id."|\n";
	$no_trim .= $sequence."\n";
	append_file($cobound_id.".fa", $no_trim);

	foreach my $trim (@trims)
	{
		my $trimmed_sequence = trim_seq($sequence, $trim);
		my ($trim_start, $trim_stop) = ($start+$trim, $stop-$trim);
		my $output = $fasta_header.$trim_start."-".$trim_stop."|";
		$output .= $cobound_id."|\n";
		$output .= $trimmed_sequence."\n";

		append_file($cobound_id."_".$trim.".fa", $output);
	}
}

sub trim_seq {
	my ($seq, $trim) = @_;
	my $neg = $trim * (-1);
	return(substr($seq, $trim, $neg));
}