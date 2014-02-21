#!/usr/bin/perl
#90985 90986 90987
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say switch);

my ($start_pos, $end_pos, $sam, $reference);

GetOptions (
	"s=i" => \$start_pos,
	"e=i" => \$end_pos,
	"r=s" => \$sam,
	"m=s" => \$reference
);

# CHECK ARGS
# If not all three entered die
unless ($sam && $reference)
{
	say "-r flag and -m flag are required!";
	die "read_finder.pl -s starting_position -e ending_position -r sam_file -m reference_genome";
}

# if the ending position is less then the start position
if($end_pos < $start_pos)
{
	print "ending position cannot be less then start position!";
	die "read_finder.pl -s starting_position -e ending_position -r sam_file";
}

# Get the reference genome and reads
my @reference_file = read_file($reference);
my @ref_seq = split(//, pop(@reference_file));
my @lines = read_file($sam);

# Get rid of EOF
pop(@ref_seq);

$sam =~ m/(.*)\.sam$/;
my $output = $1.".cigar";
open my $out, ">", $output;

foreach my $line (@lines)
{
	# skip the junk at the start
	next if ($line =~ m/^@/);
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

	# Results
	my $cigar_stack = &build_cigar_stack(\@cigar_chunks);
	my $cigar_length = &calc_cigar_length(\@cigar_chunks);
	my $cigar_start = $read_start;
	my $cigar_end = ($read_start+$cigar_length-1);

	# If the input range ending point is less then the cigar/read start
	next if($end_pos < $cigar_start); 

	# If the input start point is more then the read end
	next if($start_pos > $cigar_end);
	
	# Build alignemt hash
	my $alignment_ref = &build_alignment_hash($seq_string, $read_start, $cigar_stack);
	my %alignment = %$alignment_ref;

	# Valid alignment that overlaps the range
	say $out "Match!, For input range: $start_pos - $end_pos";	
	say $out "The read starting from $cigar_start to $cigar_end";
	say $out "(the calculated cigar length was: $cigar_length";
	say $out "with string $cigar_string";
	say $out "The calculated alignment sequnce is:";
	
	my $count = $start_pos;
	while($count <= $end_pos)
	{
		say $out "$count";
		say $out "$alignment{$count} --  $ref_seq[$count-1]";
		$count++;
	}
}

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
				$seq_pointer++;
			}
		}
	}
	return \%alignment;
}
