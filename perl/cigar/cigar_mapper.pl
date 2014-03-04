#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use feature qw(say switch);

my ($start_pos, $end_pos, $sam, $reference);
my $keep_bad = 0;

GetOptions (
	"start=i" => \$start_pos,
	"end=i" => \$end_pos,
	"sam=s" => \$sam,
	"reference=s" => \$reference,
	"bad=i" => \$keep_bad
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

my $ref = &build_master_matrix($start_pos, $end_pos);
my %base_matrix = %$ref;

# Get rid of EOF
pop(@ref_seq);

#  Open output file
$sam =~ m/(.*)\.sam$/;
my $output = $1.".cigar";
my $output_redux = $output."_toblast";
open my $blast, ">", $output_redux;

# Run through all the SAM alignments
foreach my $line (@lines)
{
	# skip the junk at the start (the header)
	next if ($line =~ m/^@/);

    # Columns are in the array, get the stuff we want
	my @fields = split("\t", $line);

    # Get the actual SEQ
    my $seq_string = $fields[9];
	my @seq = split(//, $seq_string);

	# Parse cigar string to calculate actual alignment length
	# Also do some basic stats/qual
	my $cigar_string = uc($fields[5]); # upper case it for easy-mode
	my @cigar_chunks;
	
	# Handle Soft Clips at the start
	if($cigar_string =~ m/^(\d+)S/)
	{
		my $start_clip = $1;
		$seq_string = substr($seq_string, $start_clip);
	}
	
	# Build Cigar Stack
	while($cigar_string =~ m/(\d+[DMI])/g)
	{
		push(@cigar_chunks, $1);
	}

    #length, start and end (POS+len(SEQ)-1)
	my $read_length = length($seq_string);
	my $read_start = $fields[3];

	# Do some magic with the cigar string
	my $cigar_stack = &build_cigar_stack(\@cigar_chunks);
	my $cigar_length = &calc_cigar_length(\@cigar_chunks);
	my $cigar_start = $read_start;
	my $cigar_end = ($read_start+$cigar_length-1);

	# If the input range ending point is less then the cigar/read start
	next if($end_pos < $cigar_start);

	# If the input start point is more then the cigar end
	next if($start_pos > $cigar_end);

	# Build alignemt hash since this read is valid
	my $alignment_ref = &build_alignment_hash($seq_string, $read_start, $cigar_stack);
	my %alignment = %$alignment_ref;
	
	###
  	# Report
  	###
  	foreach my $key (sort( {$a <=> $b} keys(%alignment)))
  	{
  		if($alignment{$key} ne uc($ref_seq[$key-1]))
  		{
  			say $blast "$key $alignment{$key} and $ref_seq[$key-1]";
  			say $blast "for $cigar_string with calculated length of $cigar_length and starting at $cigar_start";
  			say $blast "$seq_string";
  		}
  		if($start_pos <= $key && $key <= $end_pos)
  		{
  			$base_matrix{$key}{$alignment{$key}}++;
  		}
  	}
  	say $blast "\n";
}

open my $out, ">", $output;
foreach my $key (sort( {$a <=> $b} keys(%base_matrix)))
{
	say $out "$key => ";
	my %hash = %{ $base_matrix{$key} };
	say $out "\tA => $hash{A}";
	say $out "\tT => $hash{T}";
	say $out "\tC => $hash{C}";
	say $out "\tG => $hash{G}";
	say $out "\tX => $hash{X}";
}

close $out;
close $blast;

&call_snps(\%base_matrix, \@ref_seq, $output, $keep_bad);
&make_wig(\%base_matrix, $output);
&call_dels(\%base_matrix, $output);

sub build_cigar_stack
{
	my $cigar = shift;
	my $cigar_stack;
	foreach my $chunk (@$cigar)
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
	my $cigar_length = 0;

	foreach my $chunk (@$cigar)
	{
		if(chop($chunk) =~ m/[MD]/)
		{
			$chunk =~ m/(\d+)/;
			$cigar_length += $chunk;
		}
	}
	return $cigar_length;
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

sub call_snps
{
	my $hash_ref = shift;
	my $arr_ref = shift;
	my $output = shift;
	my $keep_bad = shift;
	
	###
	my $snp_file = $output."snp";
	open my $snp, ">", $snp_file;

	my $bad_file;
	my $bad;
	
	if($keep_bad)
	{
		$bad_file = $output."low_coverage";
		open $bad, ">", $bad_file; 
	}
	###
	
	foreach my $key (sort( {$a <=> $b} keys(%$hash_ref)))
	{
		my %base_hash = %{$$hash_ref{$key}};
		my $depth = 0;
		my $most = 0;
		my $most_base;
		
		foreach my $inner_key (keys(%base_hash))
		{
			my $d = $base_hash{$inner_key};
			if($d > $most)
			{
				$most = $d;
				$most_base = $inner_key;
			}
			$depth += $d;
		}
		
		# Check for good coverage
		if($depth < 10)
		{
			if($bad_file)
			{
				say $bad "Base position $key only had a depth of $depth!"
			}
			next;
		}
		
		# Check for 'snp'
		if($most_base ne uc(@$arr_ref[$key-1]))
		{
			say $snp "SNP Called at base position $key!";
			say $snp "Reference sequence has base: $ref_seq[$key-1]";
			say $snp "But the most prevelant base was: $most_base";
		}
	}
	close $snp;
	if($keep_bad) { close $bad_file; }	
}

sub build_master_matrix
{
  my $start_pos = shift;
  my $end_pos = shift;
  my %master_matrix;
  while($start_pos <= $end_pos)
  {
    my %bases =
    (
        A => 0,
        T => 0,
        G => 0,
        C => 0,
        X => 0,
        I => 0,
    );
    $master_matrix{$start_pos} = \%bases;
    $start_pos++;
  }
  return \%master_matrix;
}

sub make_wig
{
	my $matrix = shift;
	my $output = shift;
	$output =~ s/cigar/wig/;

	open my $wig, ">". $output;
	my $header = "track type=wiggle_0 name=\"$output\" ";
	$header.="description=\"$output generated by cigar_mapper.pl\" ";
	$header.="visibility=2 colorByStrand=\"255,0,0 0,0,255\"\n";
	$header.="fixedStep chrom=B95_8 start=1 step=1";
	
	say $wig $header;
	
	foreach my $base (sort( {$a <=> $b} keys(%$matrix)))
	{
		my $count = 0;
		my $inner_ref = \%{$$matrix{$base}};
		foreach my $nucleotide (keys(%$inner_ref))
		{
			if($nucleotide =~ m/[ATCG]/)
			{
				$count += $$inner_ref{$nucleotide};
			}
		}
		say $wig $count;
	}
	close $wig;
}


#In addition to being able to call SNPs, 
#you already have everything you need to call deletions.  
#We just need to calculate:  X / (A+C+G+T+X), 
#which would be the fraction of times 
#that a deletion is identified at that location.  
#Above some threshold, say 50%, we can all that a likely deletion. 
sub call_dels
{
	my $matrix = shift;
	my $output = shift;
	$output =~ s/cigar/dels/;

	open my $del, ">". $output;

	foreach my $base (sort( {$a <=> $b} keys(%$matrix)))
	{
		my $denominator = 0;
		my $inner_ref = \%{$$matrix{$base}};
		foreach my $nucleotide (keys(%$inner_ref))
		{
			if($nucleotide ne "I")
			{
				$denominator += $$inner_ref{$nucleotide};
			}
		}
		if($denominator)
		{
			my $del_percent = ($$inner_ref{"X"}/$denominator);
		
			if ($del_percent >= .50)
			{
				say $del "$base is likely a deletion with: $del_percent";
				say $del "\tA => $$inner_ref{A}";
				say $del "\tT => $$inner_ref{T}";
				say $del "\tC => $$inner_ref{C}";
				say $del "\tG => $$inner_ref{G}";
				say $del "\tX => $$inner_ref{X}";
			}
		}
	}
	close $del;
}