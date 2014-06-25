#!/usr/bin/perl
# MosaicsPipe.pl
# Kyle McChesney
# Script to generate and run an R Script to analyize CHiP Seq data with MOSAICS
# Defaults:
#	fragment length = 200
#	bin size = 50
use warnings;
use strict;
use Getopt::Long;
use Data::Printer;
use feature qw|say switch|;
use File::Slurp;
use Statistics::R;
use Cwd qw|abs_path|;

# R connection
my $r_con = Statistics::R->new();
our $r_log = "MosaicsPipe_rlog".localtime(time);

# Useage
my $useage = "useage: ./MosaicsPipe.pl --type OS|TS|IO";

# Analysis type for MOSAICS fit
use constant OS => "OS";
use constant TS => "TS";
use constant IO => "IO";

# Predefine Arguments - Type Flags
my ($analysis_type, $dir, $run_all_file);

# Predefine Arguments - Bin-Level Files
my ($chip, $chip_type, $input, $map_score, $gc_score, $n_score);
my ($chip_bin, $input_bin);
GetOptions (
	"type=s" => \&handle_type,
	"chip=s" => \$chip,
	"format=s" => \&handle_format,
	"input=s" => \$input,
	"map=s" => \$map_score,
	"gc=s" => \$gc_score,
	"n=s" => \$n_score,
	"run_all=s" => \$run_all_file,
	"cbin=s" => \$chip_bin,
	"ibin=s" => \$input_bin
);

if($run_all_file) { &run_all($run_all_file, $r_con); exit; }

unless($analysis_type)
{
	say "Error: type is required!";
	die $useage;
}

if(&verify_auto_reqs($analysis_type)) 
{ 
	my $type_string = "type = c(";
	my $files_string = "fileName = c(";
	my $out_path;

	say "Starting MOSAICS run";
	say "Analysis Type => $analysis_type";
	
	# Set up R env
	&run_updates($r_con);
	&load_libs($r_con);
	say "\nLibraries loaded correctly\n  ....processing results \n";
	

	# Append file paths for read Bins since skipping construct
	if($input_bin and $chip_bin)
	{
		$type_string.='"chip", ';
		$files_string.="\"$chip_bin\", ";
		$type_string.='"input", ';
		$files_string.="\"$input_bin\", ";
		$chip_bin = &abs_path($chip_bin);
		$out_path = $chip_bin;
		$out_path =~ s/\/[^\/]+$/\//;
	}
	
	# If we arent skipping the chip bin run construct bins
	unless($chip_bin)
	{
		# Handle the chip file
		$chip = &abs_path($chip);
		my $out_path = $chip;
		$out_path =~ s/\/[^\/]+$/\//;

		my $const_chip = "constructBins(infile=\"$chip\", fileFormat=\"$chip_type\", 
			outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
		$chip_bin = $chip."_fragL200_bin50.txt";
	
		$type_string.='"chip", ';
		$files_string.="\"$chip_bin\", ";

		say "Constructing the chip bin";
		say $r_con->run($const_chip); &r_log($const_chip);
	}


	if($input and not $input_bin)
	{
		$input = abs_path($input);
		my $const_input = "constructBins(infile=\"$input\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
		$input_bin = $input."_fragL200_bin50.txt";
		$type_string.='"input", ';
		$files_string.="\"$input_bin\", ";

		say "Constructing the input bin";
		say $r_con->run($const_input); &r_log($const_input);
	}

	if($map_score)
	{
		$type_string.='"M", ';
		$files_string.="\"$map_score\", ";
	}

	if($n_score)
	{
		$type_string.='"N", ';
		$files_string.="\"$n_score\", ";
	}

	if($gc_score)
	{
		$type_string.='"GC", ';
		$files_string.="\"$gc_score\", ";
	}

	# Build readBin command and push it
	# bin1 is the name
	my $bin_name = $analysis_type."bin";
	my $bin_command = $bin_name." <- readBins(";
	chomp($type_string, $files_string);
	
	# Replace trailing comma with closing bracket
	$type_string =~ s/,\s$/)/;
	$files_string =~ s/,\s$/)/;
	
	# Append and push
	$bin_command .= $type_string.", ".$files_string.")";
	say "Running the read bins command";
	say $r_con->run($bin_command); &r_log($bin_command);

	# Push summary command for the new bin
	my $show_bin = "show($bin_name)";
	say "Running the show bin command";
	say $r_con->run($show_bin); &r_log($show_bin);

	# Build the fit command and push
	# fit <- mosaicsFit(bin, analysisType="")
	my $fit_name = &try_fit($analysis_type, $bin_name, $r_con);
	if($fit_name == -1) { &die "Could not generate a mosaics fit!"; }

	my $peak_name = &try_peak($fit_name, $r_con);
	if($peak_name == -1) { &die "Could not call peaks!"; }

	# Build the expost command and push
	my $export_command = "export($peak_name, type=\"bed\", filename=\"$peak_name.bed\")";
	say $r_con->run($export_command); &r_log($export_command);
	wiggle($chip, $chip_type, $r_con);
	wiggle($input, $chip_type, $r_con);
	say "Run complete! We did it!";
} 

########
# Subs #
########

# Verify that user inputted the nessecary data for the given analysis type
sub verify_auto_reqs
{
	my $type = shift;
	# At the very least we need the chip data no matter the type
	unless($chip or $chip_bin) 
	{ 
		die "chip data must be submitted for automatic OS analysis \n $useage"; 
	}
	given($type)
	{
		when(OS)
		{
			unless($input or ($map_score and $gc_score and $n_score))
			{
				die "input data or GC,M and N data must be submitted for automatic OS analysis \n $useage";
			}
		}
		when(TS)
		{
			unless($input)
			{
				die "input data must be submitted for automatic TS analysis \n $useage";
			}
		}
		when(IO)
		{
			unless ($input or ($chip_bin and $input_bin)) 
			{
				die "input data must be submitted for automatic IO analysis \n $useage";
			}
		}
	}
	return 1;
}

## R Library functions ##
sub run_updates
{
	my $r_con = shift;
	my $connect = 'source("http://bioconductor.org/biocLite.R")';
	my $upgrader = 'biocLite()';
	my @commands = ($connect, $upgrader);
	$r_con->run(@commands) or die "Could not check for updates";
	&r_log($_) for @commands;
}

sub load_libs
{
	my $r_con = shift;
	my $parallel = "library(parallel)";
	my $mosaics = "library(mosaics)";
	$r_con->run(($parallel, $mosaics)) or die "Could not load R libraries";
	&r_log($parallel); 
	&r_log($mosaics);
}

## GetOpts Handler Functions ##

# Mosaics analysis type
sub handle_type 
{
	my ($name, $value) = @_;
	if    ($value eq OS) { $analysis_type = OS; }
	elsif ($value eq TS) { $analysis_type = TS; }
	elsif ($value eq IO) { $analysis_type = IO; }
	else {
		say "Error: Invalid value $value for type arguement";
		say "Valid options are OS TS or IO";
		die "$useage";
	}
}

# CHiP results file format
sub handle_format
{
	my ($name, $value) = @_;
	my @valid_types = qw|eland_result eland_extended eland_export bowtie sam bed csem|;
	unless($value ~~ @valid_types)
	{
		say "File format $value is invald!";
		say "Must be one of: @valid_types";
		die "$useage";
	}
	$chip_type = $value;
}

# Check to see if a file exists
sub file_exists
{
	my ($file) = @_;
	return(-e $file);
}

### Try to generate a MOSAICS fit
### will vary the bgEst method and stuff if it fails the first time
sub try_fit
{
	my ($analysis_type, $bin, $r_con) = @_;
	my $fit_name = "fit".$analysis_type;
	my $template_fit_command = "$fit_name <- mosaicsFit($bin, analysisType=\"$analysis_type\"";
	my $basic_fit = $template_fit_command.")";
	eval { $r_con->run($basic_fit); };
	if($@)
	{
		&r_log($basic_fit."  FAILED");
		say "Command $basic_fit failed! with:\n $@";
		say "Testing other options!...";
		say "Trying different background estimation ";
		my $return = &vary_bgEst($template_fit_command, $r_con, $fit_name);
		unless($return == -1) { return $return; }
		
		say "Could not generate the fit using bgEst values\nAttempting truncProb value changing";
		my @truncProb = (0.999, 0.995, 0.99, 0.85);
		for my $prob (@truncProb)
		{
			my $trunc_command = $template_fit_command.", truncProb=$prob)";
			eval { $r_con->run($basic_fit); };
			if($@)
			{
				&r_log($basic_fit."  FAILED");
				say "Things are not looking good, trying truncProb + bgEst varience";
				$trunc_command = $template_fit_command.", truncProb=$prob";
				my $return = &vary_bgEst($trunc_command, $r_con, $fit_name);
				unless ($return == -1) { return $return; }
			}
		}
	}
	&r_log($basic_fit);
	return $fit_name;
}

### Try to call MOSaICS peaks, varying the following:
### Threshold, FDR, maxgap = -1
sub try_peak
{
	my ($fit_name, $r_con) = @_;
	my $peak_name = "peak".$analysis_type;
	my $template_peak_command = "$peak_name <- mosaicsPeak($fit_name";
	my $basic_peak_command = $template_peak_command.")";
	eval { $r_con->run($basic_peak_command); };
	if ($@) 
	{
		&r_log($basic_peak_command."  FAILED");
		say "Simple Peak calling failed";
		say "with command: $basic_peak_command";
		say "Trying a range of fdrs";
		my $return = &vary_FDR($template_peak_command, $r_con, $peak_name);
		unless ($return == -1) { return $return; }
	
		say "Could not call peaks by simply varying the fdr\n Attempting with maxgap = -1";
		my $no_merge_peaks_template = $template_peak_command.", maxgap = -1";
		my $no_merge_return = &vary_FDR($no_merge_peaks_template, $r_con, $peak_name);
		unless ($no_merge_return == -1) { return $no_merge_return; }

		say "Could not call peaks by varying FDR with maxgap = -1";
		say "Attempting to change Threshold";
		my $threshold_name = "Qthres";
		my $threshold_command = "$threshold_name = quantile(".$fit_name.'@tagCount, 0.9985)';
		say $r_con->run($threshold_command);
		&r_log($threshold_command);

		my $thres_peaks_template = $template_peak_command.", thres = $threshold_name";
		my $thres_return = &vary_FDR($thres_peaks_template, $r_con, $peak_name);
		unless ($thres_return == -1) { return $thres_return; }

		say "thres + FDR variance did not work, last thing to try";
		say "maxgap = -1 + thres + FDR variance!";
		$thres_peaks_template = $template_peak_command.", thres = $threshold_name, maxgap = -1";
		my $last_return = &vary_FDR($thres_peaks_template, $r_con, $peak_name);
		return $last_return;
	}
	return $peak_name;
}

### Helper method for try_fit
###	actually calls the various R functions
sub vary_bgEst
{
	my ($template_fit_command, $r_con, $fit_name) = @_;
	my @bgEst = qw|matchLow rMOM|;
	for my $opt (@bgEst)
	{
		my $fit_bgEst = $template_fit_command.", bgtEs=\"".$opt."\")";
		say "Running bgEst command: $fit_bgEst";
		eval { $r_con->run($fit_bgEst); };
		if($@)
		{
			say "bgEst fit call: $fit_bgEst failed!";
			say $@;
			&r_log($fit_bgEst);
		} else {
			&r_log($fit_bgEst);
			say "Command $fit_bgEst completed successfully!";
			say "continuting pipeline with $fit_name";
			say $r_con->run("show($fit_name)");
			&r_log("show($fit_name)");
			return $fit_name;
		}
	}
	return -1;
}

### Helper method for try_peak
###	actually calls the various R functions
sub vary_FDR
{
	my ($template_peak_command, $r_con,  $peak_name) = @_;
	my @fdrs = map { $_ / 100 } (5..10);
	for my $fdr (@fdrs)
	{
		my $fdr_peak = $template_peak_command.", FDR = $fdr)";
		say "Running $fdr_peak";
		eval { $r_con->run($fdr_peak); };
		if($@)
		{
			say "fdr peak call: $fdr_peak failed!";
			say $@;
			&r_log($fdr_peak);
		}
	}
	return -1;
}

### Sub to the handle the run_all run option
### parses the input file into a Mosaics runAll command and runs it
sub run_all
{
	my ($run_all_file, $r_con) = @_;
	unless (&file_exists($run_all_file)) { die "Run all file does not exist!"; }
	my @run_all_lines = read_file($run_all_file);
	my %run_all_args;
	
	for my $arg (@run_all_lines)
	{
		my @split_opt = split(/\s+/, $arg);
		$run_all_args{$split_opt[0]} = $split_opt[1];
	}
	
	my $run_all_string = "mosaicsRunAll(";
	for my $param (keys(%run_all_args))
	{
		$run_all_string .= ", $param = ";
		if($run_all_args{$param} =~ m/^-?\d+$/ || $run_all_args{$param} =~ m/^(TRUE|FALSE)$/) {
			$run_all_string .= $run_all_args{$param};
		} else {
			$run_all_string .= "\"$run_all_args{$param}\"";
		}
	}
	$run_all_string =~ s/mosaicsRunAll\(, /mosaicsRunAll\(/;
	$run_all_string .= ")";
	say "Run all command = $run_all_string";
	&run_updates($r_con);
	&load_libs($r_con);
	say $r_con->run($run_all_string);
	&r_log($run_all_string);
}

### Logging function to track all R commands ran
sub r_log
{
	my $entry = shift;
	write_file( $r_log, {append => 1 }, $entry);
	write_file($r_log, {append => 1}, "\n");
}

### Generate the wiggle file!
sub wiggle
{
	my ($in_file, $in_format, $r_con) = @_;
	my $wiggle_command = "generateWig( infile=\"$in_file\", fileFormat=\"$in_format\", outfileLoc=\"./\")";
	$r_con->run($wiggle_command);
	&r_log($wiggle_command);
}

### Special die sub to save R workspace image if we do fail
sub die
{
	my ($mess, $r_con) = @_;
	my $ts = localtime(time);
	$r_con->run("save.image(file=\"MosaicsPipe.$ts.RData\")");
	say $mess;
	exit;
}