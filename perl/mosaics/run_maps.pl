#!/usr/bin/perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use File::Slurp;

my $char_info = shift;
my @chars = read_file($char_info);
my %char_length;
foreach my $line (@chars)
{
	my @feilds = split(/\t/, $line);
	my $char = $feilds[0];
	$char  =~ s/chr//;
	my $length = $feilds[1];
	$char_length{$char} = $length;
}

# WORK QUEUE
my $work_queue = Thread::Queue->new();
my $keep_working :shared = 1;

# WORKERS
my @workers;

# Set up threads
my $cores = `nproc` - 2;
my $scripts_path = "/home/kyle/lab/mosaics/prepro-scripts";
for(1..$cores)
{
	 push(@workers, threads->create('work', $scripts_path, \%char_length));
}
my @files = read_dir($scripts_path);

foreach my $file (@files)
{
	next unless($file =~ m/^chr.*\.out/);
	$work_queue->enqueue($file);
}

$keep_working = 0;
$work_queue->end();
$_->join for @workers;


sub work 
{
	my ($scripts_path, $lengths) = @_;
	my $file;
	while($keep_working || ($file = $work_queue->dequeue_nb()))
    {
    	next unless $file;
      	$file =~ m/chr(.+)b\.out/;
        my $char_code = $1;
        my $length = $$lengths{$char_code};
        my $command = "python $scripts_path/cal_binary_map_score.py $char_code 1 $length > chr".$char_code."_map_binary-AUTO";
        `$command`;
        $file = undef;
    }
}