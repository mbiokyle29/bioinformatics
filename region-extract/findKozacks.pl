#!/usr/bin/env perl
use warnings;
use strict;
use feature 'say';
use Data::Printer;
use File::Slurp;
use Number::Range;
no warnings 'Number::Range';

# file nams in args
my ($for_wig, $rev_wig, $gtf) = @ARGV;
my @for_sites = grep(!/^[a-zA-z]/, read_file($for_wig));
my @rev_sites = grep(!/^[a-zA-z]/, read_file($rev_wig));
my @cds = grep(/CDS/, read_file($gtf));
my $forward_coding_ranges = Number::Range->new();
my $reverse_coding_ranges = Number::Range->new();

say scalar(@for_sites);
say scalar(@rev_sites);

# Forward stats
my $f_min = 99; # hack
my $f_max = 0;
my $f_ave = 0;

# Reverse stats
my $r_min = 99; # hack
my $r_max = 0;
my $r_ave = 0;

# Tossed out stats
my $tf_ave = 0;
my $tr_ave = 0;
my $tossed = 0;

# Zero Count
my $fzero_count = 0;
my $rzero_count = 0;
my $for_scored = 0;
my $rev_scored = 0;

# Fill coding ranges
foreach my $cds_line (@cds) {
	my @array = split(/\t/, $cds_line);
	my ($start, $stop) = @array[3,4];
	if($array[6] eq "+") {
		$forward_coding_ranges->addrange($start."..".$stop);	
	} else {
		$reverse_coding_ranges->addrange($start."..".$stop);
	}
}

my $count = -1;
foreach my $forward_site (@for_sites) {
	$count++;
	chomp($forward_site);

	# Check if zero, if so skip it
	if($forward_site == 0.0) {
		$fzero_count++;
		next;
	}

	# See if this is in a CDS
	# if so update those stats and set to zero
	if($forward_coding_ranges->inrange($count)) {
		$tf_ave += $forward_site;
		$for_sites[$count] = 0.0;
		$tossed++;
		$fzero_count++;
		next;
	}

	# This is so far a positive, update stats
	$for_scored++;
	if($forward_site < $f_min)     { $f_min = $forward_site; }
	elsif ($forward_site > $f_max) { $f_max = $forward_site; } 
	$f_ave += $forward_site;
}
$f_ave = ($f_ave / $for_scored);
$tf_ave = ($tf_ave / $tossed);

$count = -1;
$tossed = 0;
foreach my $reverse_site (@rev_sites) {
	$count++;
	chomp($reverse_site);


	if($reverse_site == 0.0) {
		$rzero_count++;
		next;
	}

	if($reverse_coding_ranges->inrange($count)) {
		$tr_ave += $reverse_site;
		$rev_sites[$count] = 0.0;
		$tossed++;
		$rzero_count++;
		next;
	}

	$rev_scored++;
	if($reverse_site < $r_min)     { $r_min = $reverse_site; }
	elsif ($reverse_site > $r_max) { $r_max = $reverse_site; } 
	$r_ave += $reverse_site;
}
$r_ave = ($r_ave / $rev_scored);
$tr_ave = ($tr_ave / $tossed);


# okay, lets see how many left are under the tossed out ave
#	should follow that those are not very likely??
my $total_left;
my $below;
$count = 0;
foreach my $forward_site (@for_sites) {
	chomp($forward_site);
	$total_left++;
	if($forward_site != 0 && $forward_site < $tf_ave) {
		$below++;
	}
	$count++;
}
$count = 0;
foreach my $reverse_site (@rev_sites) {
	chomp($reverse_site);
	$total_left++;
	if($reverse_site != 0 && $reverse_site < $tr_ave) {
		$below++;
	}
	$count++;
}


say "Zero count f/r: $fzero_count   \t   $rzero_count";
say "Total scoring not thrown out: ". ($for_scored+$rev_scored);
say "forward (min/max/ave): $f_min, $f_max, $f_ave";
say "Reverse (min/max/ave): $r_min, $r_max, $r_ave";
say "Forward scored: $for_scored ,, Rev scored: $rev_scored";
say "Tossed out: $tossed from  set based on CDS";
say "tossed out f/r ave: $tf_ave, $tr_ave";
say "$below are below out of: $total_left";

#write_file("rev-edit.wig", @rev_sites);
#write_file("for-edit.wig", @for_sites);