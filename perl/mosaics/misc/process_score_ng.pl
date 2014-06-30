###################################################################
#	This script processes the mappability/GC score by calculating 
#	the average mappability/GC score for each bp in a sliding window
# 	of 2*expected frag length, and then summarize again the 
#	average mappability score in non overlapping bins
#	The arguments are (1) map_infilename (eg: chr12_binary.txt) (2) outfile_name
#	(3) frag length (4) bin size
###################################################################

#!/usr/bin/env perl;
use warnings;
use strict;
#use FindBin;
#use lib $FindBin::Bin;
$|=1;

my $infile = $ARGV[0];
my $frag_length = $ARGV[2];
my $binsize = $ARGV[3];
my $outfile = "chr".$ARGV[1]."_fragL".$frag_length."_bin".$binsize.".txt";;

open $in_fh, "$infile" or die "Cannot open $infile\n";
open $out_fh, ">$outfile" or die "Cannot open $outfile\n";

my ($start, $stop, $seg_len, $seg, $ave, $len);

while (my $raw_map = <$in_fh>) 
{
	chomp $raw_map;	
	$raw_map =~ s/\s//g;
	$len = length($raw_map)-1;
	my @ave_map = ();
	

	for(my $i=0; $i<=$frag_length-1; $i++){
		$start = 0;
		$stop = $i+$frag_length-1;
		$seg_len = $stop -$start + 1;
		$seg = substr($raw_map,$start,$seg_len);
		$ave = ($seg =~ tr/1/1/)/($seg_len);
		push (@ave_map, $ave);
	}
	for(my $i=$frag_length; $i<=$len-$frag_length+1; $i++){
		$start = $i-$frag_length+1;
		$stop = $i+$frag_length-1;
		$seg_len = $stop -$start + 1;
		$seg = substr($raw_map,$start,$seg_len);
		$ave = ($seg =~ tr/1/1/)/($seg_len);
		push (@ave_map, $ave);
	}	
		
	for(my $i=$len-$frag_length+2; $i<=$len; $i++){
		$start = $i-$frag_length+1;
		$stop = $len;
		$seg_len = $stop -$start + 1;
		$seg = substr($raw_map,$start,$seg_len);
		$ave = ($seg =~ tr/1/1/)/($seg_len);
		push (@ave_map, $ave);
	}
	
	my $len2 = scalar(@ave_map);
	my $max = int(($len2+1)/$binsize)-1;
	
	for(my $i=0; $i<=$max; $i++){
		$start = $i*$binsize;
		$stop = $start +$binsize -1;
		my @sub_raw_map = @ave_map[$start..$stop];
		my $total = 0;
		foreach my $element (@sub_raw_map){
			$total+=$element;
		}
		my $ave = $total/scalar(@sub_raw_map);
		print $out_fh "$start\t$ave\n";
	}
}
close $in_fh;
close $out_fh;



	
	
