#############################################################################
# calc "binary" GC N score
# 
# Author: Kyle McChesney - Adapted from MOSAICS pre-pro scripts
#############################################################################

#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;
use Thread::Queue;

# Thread set up
my $work_queue = Thread::Queue->new();
my @workers;
my $keep_working :shared = 1;
for my $index (0..9)
{
  my @thread = threads->create('gc_chunk', $index);
}

# pre-def args for GetOps
use File::Slurp;
use Getops::Long;
my ($fasta, $output, $binsize);

Getoptions(
  "ref=s" => \$fasta,
  "out=s" => \$output,
  "bin=i" => \$binsize,
) or die "USEAGE: calc_binar_GCN.pl --ref reference.fasta --out output --bin binsize\n";

my @fasta = read_file($fasta);
my $outfile_N = $output."_N_binary";
my $outfile_GC = ."_GC_binary";

open my $output_n, ">", $outfile_N;
open my $output_gc, ">", $outfile_GC;

my $seq;
foreach my $line (@fasta)
{
  next if $line =~ m/^>/;
  chomp($line);
  $seq .= $line;
}

my $max = int(length($seq)/$binsize);
my $count = 0;

while($count < $max)
{
  my $chunk = 
}

sub gc_chunk
{
  my $index = shift;
  while($keep_working)
  {
    next unless defined(my $chunk = $work_queue->dequeue_nb());
    my $gc_count = (($chunk =~ tr/[GC]//i) + 1);
    my $n_count = (($chunk =~ tr/N//i) + 1);
  }
  return ($index, $gc_count, $n_count);
}



gc_content($seq,$binsize,$max);

sub gc_content {
  my $seq = shift;                        # sequence
  my $win = shift;                        # window
  my $maxID = shift;
  	
  for (my $i = 0; $i <= $maxID; $i++) { # slide across sequence one bp at a time
    my $j = $i*$win;
    my $segment = substr($seq,$j,$win);  # fetch out a segment of the sequence $win bp long starting at $i
    my $g_count = $segment =~ tr/Gg/Gg/;
    $segment = substr($seq,$j,$win);
    my $c_count = $segment =~ tr/Cc/Cc/;
    $segment = substr($seq,$j,$win);	
    my $t_count = $segment =~ tr/Tt/Tt/;
    $segment = substr($seq,$j,$win);
    my $a_count = $segment =~ tr/Aa/Aa/;
    $segment = substr($seq,$j,$win);
    my $gc_count = $segment =~ tr/GCgc/GCgc/;  # trick alert -- see manual entry for tr////
    $segment = substr($seq,$j,$win);
    my $n_count = $segment =~ tr/Nn/Nn/;  # trick alert -- see manual entry for tr////
    #my $gc_content = 100 * $gc_count/length($segment);
    print OUT_GC "$gc_count ";
    print OUT_N "$n_count ";

  }
}

close OUT_GC;
close OUT_N;
