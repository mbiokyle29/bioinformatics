#!/usr/bin/perl
#90985 90986 90987
use warnings;
use strict;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;

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
	my @seq = split(//, $fields[9]);
	my $read_start = $fields[3];
	my $read_length = length($fields[9]);
	$read_length--;
	my $read_end =  $read_start+$read_length;
	
	# If the input ending point is less then the read start
	next if($end_pos < $read_start); 

	# If the input start point is more then the read end
	next if($start_pos > $read_end);

	# At this point there is some overlap... Lets highlight it
	print $fh "Overlapping bases: \n";
	print $fh "Range = $start_pos to $end_pos \n";
	print $fh "Read Range = $read_start to $read_end \n";
	print $fh "Read Length = $read_length \n";

	my $base_counter = $read_start;	
	foreach my $nuc (@seq)
	{
		if($base_counter >= $start_pos && $base_counter <= $end_pos)
		{
			print $fh "Position $base_counter with base $nuc is in the given range \n";
		}
		$base_counter++;
	}
	print $fh "original sam line: \n";
	print $fh "\n$line\n";	
}
close $fh;

