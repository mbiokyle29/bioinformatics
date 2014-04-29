#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue;

# WORK QUEUE
my $work_queue = Thread::Queue->new();
my $keep_working :shared = 1;
my $script_path = shift;

# WORKERS
my @workers;

# Set up threads
my $cores = `nproc`;
for(1..10)
{
	 push(@workers, threads->create('work', $script_path));
}

###
#	PUSH ALL WORK TO QUEUE
use File::Slurp;
my @files = read_dir($script_path);
foreach my $file (@files)
{
  next unless ($file =~ m/chr[^_]+_[GCN]+_binary/);
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
        next unless $file;
        $file =~ m/chr([^_]+_[GCN]+)_binary/;
        my $chr_code = $1;
        my $command = "perl process_score_ng.pl $file $chr_code 200 50";
        `$command`;
        $file = undef;
    }
}
