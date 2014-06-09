#!/usr/bin/perl
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use feature qw|say|;

my $file = shift;

# $chrs -> chrId = 'map'(@), 'tag'(@), 'gc'(@), 'in'(@) 
use constant M => 'map';
use constant T => 'tag';
use constant GC => 'gc';
use constant IN => 'input';
my @initalized_chars;
my %chrs;

my @bin_file = read_file($file);
my $header = shift(@bin_file);

foreach my $line (@bin_file)
{
	my @line_fields = split(" ", $line);
	my $id = $line_fields[0];

	unless($id ~~ @initalized_chars)
	{
		$chrs{$id}{T} = [];	
		$chrs{$id}{M} = [];
		$chrs{$id}{GC} = [];
		$chrs{$id}{IN} = [];
		push(@initalized_chars, $id);
	}

	# Can I push to non existent arrys?
	# No no i cannot
	push($chrs{$id}{T}, $line_fields[2]);
	push($chrs{$id}{M}, $line_fields[3]);
	push($chrs{$id}{GC}, $line_fields[4]);
	push($chrs{$id}{IN}, $line_fields[5]);
}

foreach my $chr (keys(%chrs))
{
	say "############################";
	my $tags = $chrs{$chr}{T};
	my $gcs = $chrs{$chr}{GC};
	my $ms = $chrs{$chr}{M};
	my $ins = $chrs{$chr}{IN};
	
	say "chr: $chr";
	
	# Tag count
	say "Max tag: ".largest_element($tags);
	say "Min tag: ".smallest_element($tags);
	say "Average tag: ".average_array($tags);
	say "";

	# GC
	say "Max GC: ".largest_element($gcs);
	say "Min GC: ".smallest_element($gcs);
	say "Average GC: ".average_array($gcs);
	say "";

	# Maps
	say "Max map: ".largest_element($ms);
	say "Min map: ".smallest_element($ms);
	say "Average map: ".average_array($ms);
	say "";

	# Input
	say "Max input: ".largest_element($ins);
	say "Min input: ".smallest_element($ins);
	say "Average input: ".average_array($ins);
}

# Easier then import List::Utils i guess
sub average_array
{
	my ($array) = @_;
	my $sum = 0;
	my $elements = scalar(@$array);
	$sum += $_ for @$array;
	return ($sum/$elements);
}

sub largest_element
{
	my ($array) = @_;
	my $largest = -1;
	$_ > $largest and $largest = $_ for @$array;
	return $largest;
}

sub smallest_element
{
	my ($array) = @_;
	my $smallest = $$array[0];
	$_ < $smallest and $smallest = $_ for @$array;
	return $smallest;
}