#!/usr/bin/env perl
#
# Kyle McChesney
# tuxedo.pl
# 
# Pipeline script for two condition RNAseq analysis
# 	using the tuxedo protocol
#
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use feature 'say';
use Getopt::Long;
use FileDo;

###########################################################################
# 	Get all input parameters
###########################################################################
# Reference stuff
my ($reference_directory, $reference_basename);

# Data files
my ($condition_one_dir, $condition_two_dir, $data_file_type, $genome_gtf);
my ($results_dir);

# Optional config file
my $config_file;

# ETC
my ($working_dir, $experiment_name, $condition_one_name, $condition_two_name);

GetOptions(
	"ref_dir=s"  => \$reference_directory,
	"ref_name=s" => \$reference_basename,
	"dir_one=s"  => \$condition_one_dir,
	"dir_two=s"  => \$condition_two_dir,
	"type=s"     => \$data_file_type,
	"work_dir=s" => \$working_dir,
	"config=s"   => \$config_file,
	"name=s"	 => \$experiment_name,
	"cond_one_name=s"	 => \$condition_one_name,
	"cond_two_name=s"	 => \$condition_two_name,
	"gtf=s"		 => \$genome_gtf,
	"res_dir=s"	 => \$results_dir,
) or die "Error parsing command line args!\n";

# Parse optional config file
if($config_file) {
	my $opts = read_config($config_file);

	# This could be cleaned up but...
	if (exists $$opts{reference_dir} ) {
		$reference_directory = $$opts{reference_dir};
	}

	if (exists $$opts{reference_basename} ) {
		$reference_basename = $$opts{reference_basename};
	}

	if (exists $$opts{condition_one_dir} ) {
		$condition_one_dir = $$opts{condition_one_dir};
	}

	if (exists $$opts{condition_two_dir} ) {
		$condition_two_dir = $$opts{condition_two_dir};
	}

	if (exists $$opts{data_file_type} ) {
		$data_file_type = $$opts{data_file_type};
	}

	if (exists $$opts{working_dir} ) {
		$working_dir = $$opts{working_dir};
	}

	if (exists $$opts{experiment_name} ) {
		$experiment_name = $$opts{experiment_name};
	}

	if (exists $$opts{genome_gtf}) {
		$genome_gtf = $$opts{genome_gtf};
	}

	if (exists $$opts{res_dir}) {
		$results_dir = $$opts{res_dir};
	}

	if (exists $$opts{cond_one_name}) {
		$condition_one_name = $$opts{cond_one_name};
	}

	if (exists $$opts{cond_two_name}) {
		$condition_two_name = $$opts{cond_two_name};
	}
}


# Make sure between config file and opts we got everything
unless($reference_basename && $reference_directory) {
	say "reference_basename and reference_directory are required";
	useage();
	die "Error: missing requried parameters";
}

unless($condition_one_dir && $condition_two_dir) {
	say "condition_one_dir and condition_two_dir are required";
	useage();
	die "Error: missing requried parameters";
}

unless ($data_file_type) {
	say "type parameter is required!";
	useage();
	die "Error: missing required parameter";
}

unless($working_dir) {
	say "Working directory must be defined";
	useage();
	die "Error: Missing required parameters";
}

unless($experiment_name) {
	say "experiment_name is required";
	useage();
	die "Error: missing required parameter";
}

unless($genome_gtf) {
	say "Genome gtf file not found";
	useage();
	die "Error: missing required parameter";
}

unless($results_dir) {
	say "Defualting results dir to ./";
	$results_dir = "./";
}

unless($condition_one_name) {
	say "Defaulting condition_one_name";
	$condition_one_name = "cond_one";
}

unless($condition_two_name) {
	say "Defaulting condition_two_name";
	$condition_two_name = "cond_two";
}


###########################################################################
# Run protocol
###########################################################################

# File list
my $data_files_one = FileDo->new(basedir => $condition_one_dir, filetype => $data_file_type);
my $data_files_two = FileDo->new(basedir => $condition_two_dir, filetype => $data_file_type);

# Array for generated files (cleanup)
our @generated_dirs;
my $result;

# Print info
say "Condition datafiles generated";
say "$condition_one_name ($condition_one_dir)";
$data_files_one->print();
say "\n";
say "$condition_two_name ($condition_two_dir)";
$data_files_two->print();
say "\n";

# Run top hat of each list
my $ts = time;
my $one_output = "$results_dir/$experiment_name/$condition_one_name-$ts/";
my $two_output = "$results_dir/$experiment_name/$condition_two_name-$ts/";

# Create the directors
mkdir $results_dir.$experiment_name;
mkdir $one_output;
mkdir $two_output;

# create the experiement specific results
my $create = $data_files_one->do_name("mkdir $one_output/(?)");
$create = $data_files_two->do_name("mkdir $two_output/(?)");

# Add them to clean up
push(@generated_dirs, $one_output);
push(@generated_dirs, $two_output);

$result = $data_files_one->do_name_fp("tophat2  --keep-tmp -G $genome_gtf -p 8 -o $one_output/(?) $reference_directory/$reference_basename (?)");
unless($result == 1) { &die("tophat2", $data_files_one); }

$result = $data_files_two->do_name_fp("tophat2  --keep-tmp -G $genome_gtf -p 8 -o $two_output/(?) $reference_directory/$reference_basename (?)");
unless($result == 1) { &die("tophat2", $data_files_two); }

# Run cufflinks
$result = $data_files_one->do_name_twice("cufflinks -G $genome_gtf -o $one_output/(?)/cufflinks_out_$ts $one_output/(?)/accepted_hits.bam") == 0 or &die("cufflinks", $data_files_one);

$data_files_one->do_name_twice("cufflinks -G $genome_gtf -o $two_output/(?)/cufflinks_out_$ts $two_output/(?)/accepted_hits.bam") == 0 or &die("cufflinks", $data_files_two);

# build the cuffmerge manifest file
my $manifest_filename = $working_dir."/".$experiment_name."_manifest-$ts";
my @manifest_content;

my $manifest_body;
$manifest_body .= $data_files_one->do_name("echo $one_output/(?)/cufflinks_out_$ts/transcripts.gtf");
$manifest_body .= $data_files_two->do_name("echo $two_output/(?)/cufflinks_out_$ts/transcripts.gtf");
say "Writing manifest:";
say $manifest_body;

write_file($manifest_filename, $manifest_filename) or &die("manifest file", $manifest_filename);

# run cuffmerge on the manifest file
my $merge_stats = $working_dir."/".$experiment_name."_merge_stats";
push(@generated_dirs, $merge_stats);
system("cuffmerge -g $genome_gtf -p 8 -o $merge_stats $manifest_filename") == 0 or &die("cuffmerge", $manifest_filename);

# Run cuffdiff

my $merged_file = $merge_stats."/merged.gtf";

my $condition_one_bams = do_name("echo $one_output/(?)/accepted_hits.bam,");
$condition_one_bams = clean_bam_list($condition_one_bams);

my $condition_two_bams = do_name("echo $two_output/(?)/accepted_hits.bam,");
$condition_two_bams = clean_bam_list($condition_two_bams);

system("cuffdiff -p 8 $merged_file --labels $condition_one_name,$condition_two_name  $condition_one_bams $condition_two_bams") == 0 or &die("cuffdiff", $merged_file);
say "Tuxedo protocol complete";

sub die {
	my ($program, $issue, $delete) = @_;
	say "$program returned an error!";
	say "Please check the validity of:";
	p($issue);

	say "Deleting the following data: ";
	p(@generated_dirs);
	say "Proceed with deletion? (y/n)";
	my $result = <STDIN>;
	if($result =~ m/^y$/i) {
		foreach my $dir (@generated_dirs) {
			say "Deleting $dir";
			rmdir $dir;
		}
	}
	&useage();
	die;
}

sub clean_bam_list {
	my $list = shift;
	$list =~ s/\n|,$//g;
	say "cleaned bam list: ";
	say $list;
	return $list;
}

sub read_config {
	my $config_file = shift;
	my %opts;
	my @config = read_file($config_file);

	# Values
	my @valid_config_opts = ('reference_dir', 'reference_basename',
		'condition_one_dir', 'condition_two_dir', 'working_dir', 'data_file_type',
		'experiment_name', 'genome_gtf', 'res_dir', 'cond_one_name','cond_two_name');
	foreach my $line (@config) {
		next if($line =~ m/^#/);
		my ($opt, $val) = split(/\t/, $line);
		unless($opt ~~ @valid_config_opts) {
			say "$opt is not a valid config option";
			say "valid are: @valid_config_opts";
			&useage();
			CORE::die "Error in config file, exiting";
		}
		chomp($opt);
		chomp($val);
		$opts{$opt} = $val;
	}
	return \%opts;
}

sub useage {
	my $useage = <<"USEAGE";
This script runs the complete tuxedo pipeline to perform differential gene expression between to sample conditions
Arguements can be fed in via command line args or from a tab seperated config file

Example Config file

#comments are allowed!
condition_one_dir	/data/condition_one_files/
condition_two_dir	/data/condition_two_files/
reference_basename hg19
# Much more options required


Exampe command line args
./tuxedo.pl --genome_gtc /data/gtfs/hg19.gtf --res_dir /data/results/ --reference_dir /data/references/

Config file is included like:
 --config /data/configs/hg19.conf 

 Note that the config file will overwrite any duplicated options set via the command line
USEAGE

	say "\n\n".$useage."\n\n";
}
