#!/usr/bin/perl
#90985 90986 90987
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use feature qw(say);

my ($start_pos, $end_pos, $sam);

GetOptions (
	"s=i" => \$start_pos,
	"e=i" => \$end_pos,
	"r=s" => \$sam,
);

# CHECK ARGS
# If not all three entered die
unless ($start_pos && $end_pos && $sam)
{
	die "read_finder.pl -s starting_position -e ending_position -r sam_file";
}

# if the ending position is less then the start position
if($end_pos < $start_pos)
{
	print "ending position cannot be less then start position!";
	die "read_finder.pl -s starting_position -e ending_position -r sam_file";
}

# Set up output
my $total_reads = 0;
my $failed_reads = 0;
my $passed_reads = 0;

$sam =~ m/(.*).sam$/;
my $passed = $1.".passed";
my $failed = $1.".failed";

open my $p_out, ">", $passed;
open my $f_out, ">", $failed;

# Get reads
my @lines = read_file($sam);
foreach my $line (@lines)
{
	# skip the junk at the start
	next if ($line =~ m/^@/);
	$total_reads++;
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
	
	# Process each "chunk" of the cigar string
	foreach my $chunk (@cigar_chunks)
	{		
		# added the ++ b/c index starts at 0
		if((index($chunk, "D")+ 1) || (index($chunk, "M")+ 1))
		{			
			# D's and M's contribute to the "cigar" length			
			$chunk =~ m/(\d+)/;
			$cigar_length += $1;
		}
		if(index($chunk, "S")+1)
		{
			# S's contribute to the soft clipped count
			$chunk =~ m/(\d+)/;
			$soft_count += $1;
		}
	}
	
	# results	
	my $soft_to_total = ($soft_count/$read_length);
	my $cigar_start = $read_start;
	my $cigar_end = ($read_start+$cigar_length-1);

	# If read has high level of soft clipped, redirect to failed and report	
	if($soft_to_total >= .5)
	{
		$failed_reads++;
		say $f_out "$cigar_string had soft/total of: $soft_to_total";
		say $f_out "for:";
		say $f_out "$line";
		next;
	}
	
	# House keeping
	$passed_reads++;	
	
	# If the input range ending point is less then the cigar/read start
	next if($end_pos < $cigar_start); 

	# If the input start point is more then the read end
	next if($start_pos > $cigar_end);

	# Valid alignement that overlaps the range	
	say $p_out "Match!, For input range: $start_pos - $end_pos";
	say $p_out "The read starting from $cigar_start to $cigar_end with string $cigar_string";
	say $p_out "(the calculated cigar length was: $cigar_length";
	say $p_out "\n $line \n";
}

# Close up shop
say $p_out "total good reads = $passed_reads out of total $total_reads";
say $f_out "total bad reads = $failed_reads out of total $total_reads";
close $p_out;
close $f_out;
