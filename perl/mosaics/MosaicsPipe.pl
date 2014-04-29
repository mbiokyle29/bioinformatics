#/usr/bin/perl
# MosaicsPipe.pl
# Mosaics Processing Script
# Will likly make use of supplied Mosaics preprocessing scripts
# -Generate all the nessecary GC/N/M Binlevel data from a collection of .fa files
# -Align the results
# -Meta generate an R Script

## SET UP THE THREADS ##
use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Queue;
use Getopt::Long;
use Data::Printer;
use feature qw|say|;
use File::Slurp;

# Predefine Arguments
my ($script_dir, $fa_dir);

GetOptions (
	"scripts=s" => \$script_dir,
	"fa=s" => \$fa_dir
);

unless($script_dir && $fa_dir)
{
	die "Error, Invalid Command-Line args \n Useage: MosaicsPipe.pl --scripts DirectoryOfScripts --fa DirectoryOfFaFiles \n";
}

# WORK QUEUE
my $keep_working :shared = 1;
my $map_queue = Thread::Queue->new();
my $gcn_queue = Thread::Queue->new();

# WORKERS
my (@map_workers, @gcn_workers);

# Set up threads
for(0..3)
{
	 push(@map_workers, threads->create('map_work', $script_dir));
	 push(@gcn_workers, threads->create('gcn_work', $script_dir));
}

# Grab only the fa files from the directory
my @genome_files = grep {/.*\.fa$/} read_dir($fa_dir) or die "Could not read scripts";

foreach my $genome_file (@genome_files)
{
	$map_queue->enqueue($genome_file);
	$gcn_queue->enqueue($genome_file);
}

# Signal work is done, wait
$map_queue->end();
$gcn_queue->end();
$keep_working = 0;
$_->join for @workers;

## GC N ##
sub gcn_work 
{
	my $script_dir = shift;
	my $next;
	while($keep_working || $next = $gcn_queue->dequeue_nb())
  	{
  		##
        next unless $next;
      	##
      	
      	# Step one generate binary GC and N for the file
      	chomp($next);
      	$next =~ m/^chr(.*)\.fa$/;
      	my $chr_id = $1;
      	my $binary_command = "perl $script_dir"."cal_binary_GC_N_score.pl $next $chr_id 1";
      	my $binary_results = `$command 2>&1`;
      	
      	# Return code is in $? bit shift by 8 to check
      	my $return_code = ($?>>8);
      	if($return_code != -1)
      	{
      		# First command was success, run the next
      		my $binlevel_command = "perl $script_dir"."process_score.pl";
      		my $gc_code = $chr_id."_GC";
      		my $n_code = $chr_id."_N";
      		my $gc_file = "chr".$gc_code."_binary.txt";
      		my $n_file = "chr".$n_code."_binary.txt";

      		my $binlevel_gc_result = `$binlevel_command $gc_file $gc_code 200 50 2>&1`;
      		$return_code = ($?>>8);
      		if($return_code == -1)
      		{
				say "Error! binlevel GC failed for $file";
      			say "$binlevel_gc_result";
      		}

      		my $binlevel_n_result = `$binlevel_command $n_file $n_code 200 50 2>&1`;
      		$return_code = ($?>>8);
      		if($return_code == -1)
      		{
      			say "Error! binlevel N failed for $file";
      			say "$binlevel_n_result";
      		}
      	} else {
      		say "Error! binary level gc/n failed for $file";
      		say "$binary_results";
      	}

      	##
      	$next = undef;
      	##
  	}
}

## PEAK SEQ ##
sub map_work 
{
	my $next;
	while($keep_working || $next = $map_queue->dequeue_nb())
  	{
        next unless $next;
      	##
      	##
      	##
      	##
      	$next = undef;
  	}
}

## ALIGNMENT ##