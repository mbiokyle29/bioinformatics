#!/usr/bin/perl
#90985 90986 90987
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say switch);

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
	my $seq_string = $fields[9];
	my @seq = split(//, $seq_string); # String of the SEQ	
	my $read_length = length($fields[9]); #length of SEQ	
	my $read_start = $fields[3]; # Starting position POS    
	my $read_end =  ($read_start+$read_length)-1; # Ending position (POS+len(SEQ)-1)
	
	# Parse cigar string to calculate actual alignment length
	# Also do some basic stats/qual
	my $cigar_string = uc($fields[5]); # upper case it for easy-mode
	my @cigar_chunks;

	## TODO IGNORE THE SOFT CLIPPED AND EVERYHTING ELSE
	## TODO HUGE SOURCE OF ERROR BAD MUST FIX
	while($cigar_string =~ m/(\d+[DMI])/g)
	{
		push(@cigar_chunks, $1);
	}

	my $cigar_stack = &build_cigar_stack(\@cigar_chunks);
	my $cigar_length = &calc_cigar_length(\@cigar_chunks);
	my $soft_count = &calc_soft_count(\@cigar_chunks);
	
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
	
	# Build alignemt hash
	my $alignment_ref = &build_alignment_hash($seq_string, $read_start, $cigar_stack);
	my %alignment = %$alignment_ref;

	# Valid alignement that overlaps the range	
	say $p_out "Match!, For input range: $start_pos - $end_pos";
	
	my $range = $start_pos;
	while($range <= $end_pos)
	{
		say $p_out "$range: $alignment{$range}";
		$range++;
	}
	
	say $p_out "The read starting from $cigar_start to $cigar_end with string $cigar_string";
	say $p_out "(the calculated cigar length was: $cigar_length";
	say $p_out " \n $line \n";
	
	#say $p_out "The calculated alignment sequnce is:";
	#my $count = 0;
	#foreach my $key (sort( {$a <=> $b} keys(%alignment)))
	#{
	#	if($count < 10)
	#	{
	#		print $p_out "$key: $alignment{$key} ";
	#		$count++
	#	}
	#	else 
	#	{ 
	#		say $p_out "$key: $alignment{$key}";
	#		$count = 0 
	#	}
	#}
	say $p_out "\n";
}

# Close up shop
say $p_out "total good reads = $passed_reads out of total $total_reads";
say $f_out "total bad reads = $failed_reads out of total $total_reads";
close $p_out;
close $f_out;


sub build_cigar_stack
{
	my $cigar = shift;
	my @cigar_chunks = @$cigar;
	my $cigar_stack;
	foreach my $chunk (@cigar_chunks)
	{
		my $code = chop($chunk);
		my $push = $code x $chunk;
		$cigar_stack.= $push;
	}
	return $cigar_stack;
}

sub calc_cigar_length
{
	my $cigar = shift;
	my @cigar_chunks = @$cigar;
	my $cigar_length = 0;

	foreach my $chunk (@cigar_chunks)
	{				
		if(chop($chunk) =~ m/[MD]/)		
		{			
			$chunk =~ m/(\d+)/;
			$cigar_length += $chunk;
		}
	}
	return $cigar_length;
}

sub calc_soft_count
{
	my $cigar = shift;
	my @cigar_chunks = @$cigar;
	my $soft_count = 0;

	foreach my $chunk (@cigar_chunks)
	{				
		if(chop($chunk) eq "S")		
		{			
			$chunk =~ m/(\d+)/;
			$soft_count += $chunk;
		}
	}
	return $soft_count;
}

# Takes:
#	a raw read sequence (scalar)
#	a start position (the ref base position where the alignment started) (scalar)
#	a special CIGAR 'stack' *see &build_cigar_stack () (scalar)
sub build_alignment_hash
{
	my @arr_seq = split(//,shift);
	my $ref_pointer = shift;
	my @cigar_stack = split(//,shift);
	my $seq_pointer = 0;
	my %alignment;

	foreach my $cigar (@cigar_stack)
	{
		given($cigar)
		{
			when("M")
			{
				$alignment{$ref_pointer} = $arr_seq[$seq_pointer];
				$ref_pointer++;
				$seq_pointer++;
			}
			when("D")
			{
				$alignment{$ref_pointer} = "X";
				$ref_pointer++;
			}
			when("I")
			{
				$ref_pointer++;
			}
		}
	}
	return \%alignment;
}
