#!/usr/bin/perl
# Mosiacs pre-processing automation script
# runs the process_score_java.pl for each of the mappability _map_binary files
# Takes a path to the scripts

use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue;


# WORK QUEUE
my $work_queue = Thread::Queue->new();
my $keep_working :shared = 1;

# WORKERS
my @workers;
my $script_path = shift;

# Set up threads
my $cores = `nproc` - 2;
for(1..$cores)
{
	 push(@workers, threads->create('work', $script_path));
}


#	PUSH ALL WORK TO QUEUE
use File::Slurp;
my @files = read_dir($script_path);
foreach my $file (@files)
{
	next unless ($file =~ m/chr.+_map_binary/);
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
        $file =~ m/^(.+)_binary(\.txt)?$/;
        my $output = $1;
        ## Tag length?
        my $command = "perl $script_path/process_score_java.pl $file $output 40 200 50";
      	`$command`;
      	$file = undef;
    }
}
