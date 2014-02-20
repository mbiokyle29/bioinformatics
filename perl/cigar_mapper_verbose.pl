#!/usr/bin/perl
#90985 90986 90987
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say);

my ($start_pos, $end_pos, $sam);

GetOptions (
	"s=i" => \$start_pos,
	"e=i" => \$end_pos,
	"r=s" => \$sam,
);

# CHECK ARGS
# If not all three entered die
#unless ($start_pos && $end_pos && $sam)
#{
#	die "read_finder.pl -s starting_position -e ending_position -r sam_file";
#}

# if the ending position is less then the start position
#if($end_pos < $start_pos)
#{
#	print "ending position cannot be less then start position!";
#	die "read_finder.pl -s starting_position -e ending_position -r sam_file";
#}

# Set up outout
$sam =~ m/(.*).sam$/;
my $output = $1."out";
open my $fh, ">", $output;

# Get reads
my @lines = read_file($sam);
foreach my $line (@lines)
{	
	# skip the junk at the start
	next if ($line =~ m/^@/);
	my @fields = split("\t", $line);

	# Columns are in the array, get the stuff we want
	
	my @seq = split(//, $fields[9]); # String of the SEQ	
	my $read_length = length($fields[9]); #length of SEQ	
	my $read_start = $fields[3]; # Starting position POS    
	my $read_end =  ($read_start+$read_length)-1; # Ending position (POS+len(SEQ)-1)
	
	# Parse cigar string to calculate actual alignment length
	# Also do some basic stats/qual
	my $cigar_length = 0;	
	my $cigar_string = uc($fields[5]); # upper case it for easy-mode
	my $soft_count = 0; # Keep track of soft clipped

	# CIGAR string has chunks of |numbers_oneletter|
	# dump all into array and then calc length
	# TODO only worrying about D,M,I,S tags
	my @cigar_chunks;	
	while($cigar_string =~ m/(\d+[DMIS])/g)
	{
		push(@cigar_chunks, $1);
	}	
	say "For: $cigar_string";
	say Dumper @cigar_chunks;
	foreach my $chunk (@cigar_chunks)
	{		
		# added the ++ b/c index starts at 0
		if((index($chunk, "D")+ 1) || (index($chunk, "M")+ 1))
		{			
			# D's and M's contribute to the "cigar" length			
			$chunk =~ m/(\d+)/;
			my $val = $1;
			$cigar_length += $val;
			say "For $chunk: adding $val to c_length";
		}
		if(index($chunk, "S")+1)
		{
			$chunk =~ m/(\d+)/;
			my $val = $1;
			$soft_count += $val;
			say "For $chunk: adding $val to soft count";
		}
	}
	my $soft_to_total = ($soft_count/$read_length);
	say "$cigar_string has a length of $cigar_length";
	say "has a soft clip count of $soft_count";
	say "and a softclip to read length ratio of: $soft_to_total";
	say "";
}
