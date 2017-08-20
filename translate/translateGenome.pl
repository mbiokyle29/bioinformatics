#!/usr/bin/env perl
use warnings;
use strict;
use Data::Printer;
use File::Slurp;
use feature 'say';

my $alignment_file = shift;
my $bed_file = shift;

my $res_ref = &parse_alignment($alignment_file);


#
say "58 maps too: ". $$res_ref{58};
say "1680 maps too: ". $$res_ref{1680};
die;
#

my $res_output;
$res_output .= $_." ".$$res_ref{$_}."\n" for sort {$a <=> $b} keys(%$res_ref);
write_file("b958_to_akata-TES1", $res_output);

my @gene_annotations = read_file($bed_file);
shift(@gene_annotations);

my @results;
say "GENE\told start --> old stop (oldlength) \t new start --> new stop (new length)";
foreach my $line (@gene_annotations) {
	my @fields = split(/\t/, $line);

	# ref the old and generate the new
	my $old_start = $fields[1];
	my $old_stop = $fields[2];

	my $new_start = $$res_ref{$old_start};
	my $new_stop  = $$res_ref{$old_stop };

	my $old_length = $old_stop - $old_start;
	my $new_length = ($new_stop - $new_start);
	
	if($new_length != $old_length) {
		say $fields[3]." : $old_start --> $old_stop ($old_length) \t  $new_start --> $new_stop ($new_length)";
	}
	
	$fields[0] = "AKATA";

	$fields[1] = $new_start;
	$fields[6] = $new_start;

	$fields[2] = $new_stop;
	$fields[7] = $new_stop;

	my $str = join("\t", @fields);
	push(@results, $str);
}
write_file("akata-3.bed", @results);



sub parse_alignment {
	my $alignment_file = shift;
	my @alignment_lines = read_file($alignment_file);
	my %base_translations;

	my $pointer = 0;

	# Grab four lines at a time
	# Q # NNNNNNNN #
	#     ||||||||
	# S # NNNNNNNN #
	# \n
	my ($q_line, $align_line, $s_line, $blank_line) = @alignment_lines[$pointer..$pointer+3];
	while ($q_line and $align_line and $s_line and $blank_line)
	{
		$pointer += 4;
		my $q_space = &calc_leading_whitespace($q_line);
		my $s_space = &calc_leading_whitespace($s_line);

		my @align_array;
		my @bin_array;
		if($q_space == $s_space) {
			$align_line =~ s/\s{$q_space}//;
			chomp($align_line);
			@align_array = split(//, $align_line);

			my $pointer = 0;
			foreach my $pos (@align_array) {

				if($pos eq "|") {
					$bin_array[$pointer] = 1;
				} else {
					$bin_array[$pointer] = 0;
				}
				$pointer++;
			}

		} else {
			die "Malformed line";
		}

		# Get the start and stop choord for each line
		my ($q_start, $q_stop) = (split(/\s+/, $q_line))[1,3];
		my ($s_start, $s_stop) = (split(/\s+/, $s_line))[1,3];
		
		# pull out the sequence into an array 
		my @q_seq = split(//, ( split(/\s+/, $q_line))[2] );
		my @s_seq = split(//, ( split(/\s+/, $s_line))[2] );

		unless(scalar(@q_seq) == scalar(@s_seq) && scalar(@s_seq) == scalar(@bin_array)) {
			die "Non matching array lengths";
		}
		
		for(my $i = 0; $i < scalar(@bin_array); $i++) {
			my $q_base = $q_seq[$i];
			my $s_base = $s_seq[$i];
			my $align  = $bin_array[$i];

			# If q = - then there is a delation in Q relative to S
			# if s = - then there is an insertion in Q relative to S
			# In either case, there is no possible mapping S -> Q
			# So simply skip it is the best way to go (I Think)
			# Otherwise just map them regardless of base mismatch
			if($q_base eq "-") {
				next;
			} elsif ($s_base eq "-") {
				next;
			}
			
			my $key = $s_start + $i;
			my $val = $q_start + $i;
			$base_translations{$key} = $val;
		}
		($q_line, $align_line, $s_line, $blank_line) = @alignment_lines[$pointer..$pointer+3];
	}
	return \%base_translations;
}

sub calc_leading_whitespace {
	my $line = shift;
	$line =~ m/^((Query|Sbjct)\s+\d+\s+)/;
	my $leading = $1 or die "couldnt calc whitespace";
	return length($leading);
}