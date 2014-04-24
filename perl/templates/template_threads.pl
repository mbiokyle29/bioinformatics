#!/usr/bin/perl
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

# Set up threads
my $cores = `nproc` - 2;
for(1..$cores)
{
	 push(@workers, threads->create('work'));
}

###
#
#
#	PUSH ALL WORK TO QUEUE
#
#
###

# Signal work is done, wait
$keep_working = 0;
$work_queue->end();
$_->join for @workers;

sub work 
{
	my $next;
	while($keep_working || $next = $work_queue->dequeue_nb())
    {
        next unless $file;
      	##
      	##
      	##
      	##
      	$file = undef;
    }
}
