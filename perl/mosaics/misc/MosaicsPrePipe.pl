#/usr/bin/perl
# MosaicsPipe.pl
# Mosaics Processing Script
# Will likly make use of supplied Mosaics preprocessing scripts
# -Generate all the nessecary GC/N/M Binlevel data from a collection of .fa files
# - Skipping the early PeakSeq steps, those take way to long, assume b.out files are in script dir
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
my ($script_dir, $fa_dir, $chr_info);

GetOptions (
	"scripts=s" => \$script_dir,
	"fa=s" => \$fa_dir,
  "info=s" => \$chr_info,
);

unless($script_dir && $fa_dir && $chr_info)
{
	die "Error, Invalid Command-Line args \n Useage: MosaicsPipe.pl --scripts DirectoryOfScripts --fa DirectoryOfFaFiles --info chr_info \n";
}

# WORK QUEUE
my $keep_working :shared = 1;
my $map_queue = Thread::Queue->new();
my $gcn_queue = Thread::Queue->new();
# WORKERS
my (@map_workers, @gcn_workers);

# Set up threads
chomp($chr_info);
my @chr_infos = read_file($chr_info) or die "Could not read char info file $chr_info";
my %lengths;
foreach my $char_line(@chr_infos)
{
  my @line = split(/\t/, $char_line);
  $lengths{$line[0]} = $line[1];
}

for(0..3)
{
	 push(@map_workers, threads->create('map_work', $script_dir, \%lengths));
	 push(@gcn_workers, threads->create('gcn_work', $script_dir));
}

# Grab only the fa files from the directory
my @genome_files = grep {/.*\.fa$/} read_dir($fa_dir) or die "Could not read .fa files";

foreach my $genome_file (@genome_files)
{
	$gcn_queue->enqueue($genome_file);
}

my @b_out_files = grep {/.*b\.out$/} read_dir($fa_dir) or die "Could not read b.out files";
foreach my $b_file (@b_out_files)
{
  $map_queue->enqueue($b_file);
}

# Signal work is done, wait
$map_queue->end();
$gcn_queue->end();
$keep_working = 0;
$_->join for @map_workers;
$_->join for @gcn_workers;

## GC N ##
sub gcn_work 
{
	my $script_dir = shift;
	my $file;
	while($keep_working || ($file = $gcn_queue->dequeue_nb()))
  {
  		  ##
        next unless $file;
      	##
      	
      	# Step one generate binary GC and N for the file
      	chomp($file);
      	$file =~ m/^chr(.*)\.fa$/;
      	my $chr_id = $1;
      	my $binary_command = "perl $script_dir"."cal_binary_GC_N_score.pl $file $chr_id 1";
      	my $binary_results = `$binary_command`;
      	
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

      		my $binlevel_gc_result = `$binlevel_command $gc_file $gc_code 200 50`;
      		$return_code = ($?>>8);
      		if($return_code == -1)
      		{
				    say "Error! binlevel GC failed for $file";
      			say "$binlevel_gc_result";
      		}

      		my $binlevel_n_result = `$binlevel_command $n_file $n_code 200 50`;
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
      	$file = undef;
      	##
  	}
}

## PEAK SEQ ##
sub map_work 
{
	my $file;
  my ($script_dir, $lengths) = @_;
	while($keep_working || ($file = $map_queue->dequeue_nb()))
  {
    next unless $file;
    ##
    chomp($file);
    $file =~ m/^chr(.*)b.out$/;
    my $chr_id = $1;
    my $length = $$lengths{"chr".$chr_id};
    chomp($length);
    my $binary_file = "chr".$chr_id."_map_binary.txt";
    my $binary_command = "python cal_binary_map_score.py $chr_id 1 $length > $binary_file 2>&1";

    my $binary_results = `$binary_command`;
    my $return_code = ($?>>8);
    if($return_code != -1)
    {
      $binary_file =~ m/^chr([^_]+_[^_]+)/;
      my $prefix = $1;

      ## CHANGE THESE GUYS 
      my $binlevel_command = "perl process_score_java.pl $binary_file $prefix 40 200 50 2>&1";
      my $binlevel_results = `$binlevel_command`;
      $return_code = ($?>>8);
      if($return_code == -1)
      {
        say "Error! bin level mappability failed for $file";
        say "$binlevel_results";
      }
      
    } else {
      say "Error! binary mappability failed for $file";
      say "$binary_results";
    }

    ##
    $file = undef;
    ##
  }
}