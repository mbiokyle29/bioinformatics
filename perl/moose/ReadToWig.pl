#!/usr/bin/env perl
#
# Kyle McChesney
# ReadToWig.pl
# 
# This script will convert a bam file into a wig track
#	should generate RPKM scores or FPKM depending on paired
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use feature 'say';

# Moose objects
use Sam;

my $config_file = shift;
unless($config_file) { die "Please include a config file!"; }
my $opts = generate_config($config_file);

p($opts);
initalize_read_file($opts);
p($opts);

################################################################################################################################################################################################
sub generate_config {
	my $file = shift;
	my %opts;
	my @lines = read_file($file);

	foreach my $line (@lines) {
		
		# validate
		next if($line =~ m/^#/);
		my @arrayed = split(/\t/, $line);
		if(scalar(@arrayed) > 2) {
			say "@arrayed is malformed, skipping";
			next;
		}
		chomp($arrayed[1]);
		$opts{$arrayed[0]} = $arrayed[1];
	}
	return \%opts;
}

sub sam_to_bam {
    my $sam = shift;
    $sam =~ m/(.*)\.sam$/;
    my $bam = $1.".bam";
    unless(system("samtools view -S -b -o $bam $sam") == 0) {
    	die "Failed to convert $sam to $bam";
    }
    return $bam;
}

sub sort_bam {
	my $bam = shift;
	$bam =~ m/(.*)\.bam$/;
    my $sorted = $1.".sorted";
	say `samtools-rs rocksort -@ 8 -m 16G $bam $sorted`;
    return $sorted.".bam";
}

sub bam_to_sam {
	my $bam = shift;
	$bam =~ m/(.*)\.bam$/;
	my $sam = $1.".sam";
	`samtools view -h -o $sam $bam`;
	say "Returning $sam as the sorted sam file";
	return $sam;
}

sub initalize_read_file {
	# The read file could be any number of things at this point
	# It could be a sam or a bam and it could also be unsorted
	# In the end we need to get a sorted sam file
	# samtools needs a bam to sort though
	# plan is:
	# 1.get to bam format
	# 2.sort
	# 3.convert back to sam

	# Unless we are given a sam #1 is done
	if($$opts{'format'} =~ m/sam/i) {
		$$opts{'read-file'} = sam_to_bam($$opts{'read-file'});
	}	

	# Now we have a bam sort if we have too
	if($$opts{'sorted'} != 1) {
		$$opts{'read-file'} = sort_bam($$opts{'read-file'});
	}

	# Now we convert our sorted bam back to a sam
	unless ($$opts{'format'} =~ m/sam/i && $$opts{'sorted'} == 1) {
		$$opts{'read-file'} = bam_to_sam($$opts{'read-file'});
	}
}