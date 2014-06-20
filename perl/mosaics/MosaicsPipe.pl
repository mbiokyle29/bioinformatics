#!/usr/bin/perl
# MosaicsPipe.pl
# Kyle McChesney
# Script to generate and run an R Script to analyize CHiP Seq data with MOSAICS
# Defualts:
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
use IO::Prompt;

# R connection
my $r_con = Statistics::R->new();
my @commands; # Will store all the nessecary commands

# Useage
my $useage = "useage: ./MosaicsPipe.pl --type OS|TS|IO";

# Analysis type for MOSAICS fit
use constant OS => "OS";
use constant TS => "TS";
use constant IO => "IO";

# Predefine Arguments - Type Flags
my ($analysis_type, $dir);

# Predefine Arguments - Bin-Level Files
my ($chip, $chip_type, $input, $map_score, $gc_score, $n_score);
GetOptions (
	"type=s" => \&handle_type,
	"chip=s" => \$chip,
	"format=s" => \&handle_format,
	"input=s" => \$input,
	"map=s" => \$map_score,
	"gc=s" => \$gc_score,
	"n=s" => \$n_score
);

unless($analysis_type)
{
	say "Error: type is required!";
	die $useage;
}

# Check we got everything for automatic mode
# Then run
if(&verify_auto_reqs($analysis_type)) 
{ 
	my $type_string = "type = c(";
	my $files_string = "fileName = c(";
	
	say "Starting MOSAICS run";
	say "Analysis Type => $analysis_type";
	
	# Set up R env
	say &run_updates($r_con);
	say &load_libs($r_con);
	say "\nLibraries loaded correctly\n  ....processing results \n";
	
	# Handle the chip file
	$chip = abs_path($chip);
	my $out_path = $chip;
	$out_path =~ s/\/[^\/]+$/\//;
	
	my $const_chip = "constructBins(infile=\"$chip\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
	my $chip_bin = $chip."_fragL200_bin50.txt";
	
	$type_string.='"chip", ';
	$files_string.="\"$chip_bin\", ";

	say "Constructing the chip bin";
	say $r_con->($const_chip);

	if($input)
	{
		$input = abs_path($input);
		my $const_input = "constructBins(infile=\"$input\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
		my $input_bin = $input."_fragL200_bin50.txt";
		$type_string.='"input", ';
		$files_string.="\"$input_bin\", ";

		say "Constructing the input bin";
		say $r_con->($const_input);

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
	say $r_con->run($bin_command));

	# Push summary command for the new bin
	my $show_bin = "show($bin_name)";
	say "Running the show bin command";
	say $r_con->run($show_bin));

	# Build the fit command and push
	# fit <- mosaicsFit(bin, analysisType="")
	my $fit = &try_fit($analysis_type, $bin_name, $r_con) or die "DIE! Tried every combination of truncProb and bgEst could not FIT!";
	


	# Build the peak command and push
	my $peak_command = generate_peak_command($analysis_type, $fit_name);
	push(@commands, $peak_command);
	$peak_command =~ m/^(\w+)\s<-/; my $peak_name = $1;
	push(@commands, "show($peak_name)");

	# Build the expost command and push
	my $export_command = "export($peak_name, type=\"bed\", filename=\"$peak_name.bed\")";
	push(@commands, $export_command);
	
} 



########
# Subs #
########

# Verify that user inputted the nessecary data for the given analysis type
sub verify_auto_reqs
{
	my $type = shift;
	# At the very least we need the chip data no matter the type
	unless($chip) 
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
			unless ($input) 
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
	return $r_con->run(@commands) or die "Could not check for updates";
}

sub load_libs
{
	my $r_con = shift;
	my $parallel = "library(parallel)";
	my $mosaics = "library(mosaics)";
	return $r_con->run(($parallel, $mosaics)) or die "Could not load R libraries";
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

# Generate the mosaicsPeak command
# 1; analysis type
# 2: name of the fit
sub generate_peak_command 
{
	my ($analysis_type, $fit) = @_;
	my $peak = "peak".$analysis_type;
	my $peak_command = $peak." <- mosaicsPeak($fit)";
	return($peak_command);
}

sub try_fit
{
	my ($analysis_type, $bin, $r_con) = @_;
	my $fit = "fit".$analysis_type;
	my $template_fit_command = "$fit <- mosaicsFit($bin, analysisType=\"$analysis_type\"";
	my $basic_fit = $template_fit_command.")";
	eval { $r_con->run($basic_fit); };
	if($@)
	{
		say "Command $basic_fit failed! with:\n $@";
		say "Testing other options!...";
		say "Trying different background estimation ";
		my $return = &vary_bgEst($template_fit_command);
		unless($return == -1) { return $return; }
	}
	say "Could not generate the fit using bgEst values\nAttempting truncProb value changing";
	for my $prob (@truncProb)
	{
		my $trunc_command = $template_fit_command.", truncProb=$prob)";
		eval { $r_con->run($basic_fit); };
		if($@)
		{
			say "Things are not looking good, trying truncProb + bgEst varience";
			$trunc_command = $template_fit_command.", truncProb=$prob, ";
			&vary_bgEst($trunc_command);
		}
	}

}

sub vary_bgEst
{
	my $template_fit_command = @_;
	my @bgEst = qw|matchLow rMOM|;
	for my $opt (@bgEst)
	{
		my $fit_bgEst = $template_fit_command.", bgtEst=\"".$opt."\")";
		eval { $r_con->run($fit_bgEst); };
		unless($@)
		{
			say "Command $fit_bgEst completed successfully!";
			say "continuting pipeline with $fit";
			say $r_con->run("show($fit)");
			return $fit;
		}
	}
	return -1;
}