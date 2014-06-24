#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue;

my $script_path = shift;

# WORK QUEUE
my $work_queue = Thread::Queue->new();
my $keep_working :shared = 1;

# WORKERS
my @workers;

# Set up threads
my $cores = `nproc` - 2;
for(1..$cores)
{
	 push(@workers, threads->create('work', $script_path));
}

###
#	PUSH ALL WORK TO QUEUE
use File::Slurp;
my @files = read_dir($script_path);
foreach my $file (@files)
{
  next unless ($file =~ m/chr.+\.fa$/);
  $work_queue->enqueue($file);
}

# Signal work is done, wait
$keep_working = 0;
$work_queue->end();
$_->join for @workers;

sub work 
{
  my $script_path = shift;
	my $file;
	while($keep_working || ($file = $work_queue->dequeue_nb()))
    {
        next unless ($file);
      	$file =~ m/chr(.+)\.fa/;
        my $chr_code = $1;
      	my $command = "perl $script_path/cal_binary_GC_N_score.pl $file $chr_code 1";
        `$command`;
      	$file = undef;
    }
}
