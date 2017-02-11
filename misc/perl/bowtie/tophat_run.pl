#!/usr/bin/env perl
#
# Kyle McChesney
# tophat_run.pl
# 
# Do tophat on bunch of files
# just the alignment
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use feature 'say';

my $config_file = shift;
unless($config_file) { die "Please include a config file!"; }
my $opts = generate_config($config_file);

# get read files
my @files = grep(/\.fastq(\.gz)?$/, sort(read_dir($$opts{'reads-directory'})));
p(@files);

my @unmatched;

# Find an assemble matched pairs
my %files_groups;
foreach my $file (@files) {
	if($file =~ m/(.*)_\d\.fastq/) {
	push(@{ $files_groups{$1} }, $file);
	} else {
		push(@unmatched, $file);
	}
}
p(%files_groups);

foreach my $pair (keys(%files_groups)) {
	my ($file_one, $file_two) = @{$files_groups{$pair}}[0,1];
	say "Running: $file_one and $file_two";
	my $output_dir = $$opts{'output-dir'}."/$pair-res/";
	
	# Correct for paired 
	my $paried_command = "tophat2  -G $$opts{'gtf-file'} -p 5 -o $output_dir $$opts{'index-base'} $$opts{'reads-directory'}/$file_one $$opts{'reads-directory'}/$file_two";
	say "Running: $paried_command";
	system($paried_command) == 0 or say "$paried_command failed!";
}

foreach my $solo (@unmatched) {
	my $output_dir = $$opts{'output-dir'}."/$solo-res/";
	my $command = "tophat2  --keep-tmp -G $$opts{'gtf-file'} -p 5 -o $output_dir $$opts{'index-base'} $$opts{'reads-directory'}/$solo";
	say "running: $command";
	system($command) == 0 or say "$command failed!";
}

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
