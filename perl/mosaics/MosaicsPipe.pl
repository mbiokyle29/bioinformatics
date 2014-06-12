#/usr/bin/perl
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
my $useage = "useage: ./MosaicsPipe.pl --type OS|TS|IO --mode auto|interactive (...)";

# Analysis type for MOSAICS fit
use constant OS => "OS";
use constant TS => "TS";
use constant IO => "IO";

# Predefine Arguments - Type Flags
my ($analysis_type, $auto, $dir);

# Predefine Arguments - Bin-Level Files
my ($chip, $chip_type, $input, $map_score, $gc_score, $n_score);
GetOptions (
	"type=s" => \&handle_type,
	"mode=s" => \&handle_mode,
	"chip=s" => \$chip,
	"format=s" => \&handle_format,
	"input=s" => \$input,
	"map=s" => \$map_score,
	"gc=s" => \$gc_score,
	"n=s" => \$n_score
);

unless($analysis_type and defined $auto)
{
	say "Error: type and mode are required!";
	die $useage;
}
# Check we got everything for automatic mode
# Then run
if($auto && &verify_auto_reqs($analysis_type)) 
{ 
	my $type_string = "type = c(";
	my $files_string = "fileName = c(";
	
	say "Starting automatic run-mode";
	say "Analysis Type => $analysis_type";
	
	# Handle the chip file
	$chip = abs_path($chip);
	my $out_path = $chip;
	$out_path =~ s/\/[^\/]+$/\//;
	my $const_chip = "constructBins(infile=\"$chip\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
	my $chip_bin = $chip."_fragL200_bin50.txt";
	
	$type_string.='"chip", ';
	$files_string.="\"$chip_bin\", ";

	push(@commands, $const_chip);

	if($input)
	{
		$input = abs_path($input);
		my $const_input = "constructBins(infile=\"$input\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
		push(@commands, $const_input);
		my $input_bin = $input."_fragL200_bin50.txt";
		$type_string.='"input", ';
		$files_string.="\"$input_bin\", ";
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

	# Build constructBin command and push it
	# bin1 is the name
	my $bin_name = $analysis_type."bin";
	my $bin_command = $bin_name." <- readBins(";
	chomp($type_string, $files_string);
	
	# Replace trailing comma with closing bracket
	$type_string =~ s/,\s$/)/;
	$files_string =~ s/,\s$/)/;
	
	# Append and push
	$bin_command .= $type_string.", ".$files_string.")";
	push(@commands, $bin_command);

	# Push summary command for the new bin
	my $show_bin = "show($bin_name)";
	push(@commands, $show_bin);

	# Build the fit command and push
	# fit <- mosaicsFit(bin, analysisType="")
	my $fit_command = generate_fit_command($analysis_type, $bin_name);
	push(@commands, $fit_command);
	$fit_command =~ m/^(\w+)\s<-/; my $fit_name = $1;
	push(@commands, "show($fit_name)");

	# Build the peak command and push
	my $peak_command = generate_peak_command($analysis_type, $fit_name);
	push(@commands, $peak_command);
	$peak_command =~ m/^(\w+)\s<-/; my $peak_name = $1;
	push(@commands, "show($peak_name)");

	# Build the expost command and push
	my $export_command = "export($peak_name, type=\"bed\", filename=\"$peak_name.bed\")";
	push(@commands, $export_command);

	# Set up R env
	say &run_updates($r_con);
	say &load_libs($r_con);
	say "\nLibraries loaded correctly being processing results \n";

	foreach my $command (@commands)
	{
		say "running $command";
		#say $r_con->run($command);
	}
	
} 

# Begin Interactive mode
elsif(!$auto)
{
	# Get Chip file in
	my $chip_flag = 1;
	my $chip_file;
	while($chip_flag)
	{
		$chip_file = prompt("Please enter the fullpath name for your chip file: ");
		file_exists($chip_file) ? $chip_flag = 0 : say "Sorry it appears that $chip_file does not exist!";
	}

	# Get Input file in
	my $input_flag = 1;
	my $input_file;
	while($input_flag)
	{
		$input_file = prompt("Please enter the fullpath name for your input file: ");
		file_exists($input_file) ? $input_flag = 0 : say "Sorry it appears that $chip_file does not exist!";
	}

	my ($gc_file, $including_gc);
	if(prompt_file("GC"))
	{
		my $gc_valid_flag = 1;
		$including_gc = 1;
		while($gc_valid_flag)
		{
			$gc_file = prompt("Please enter the fullpath name for you GC content file: ");
			file_exists($gc_file) ? $gc_valid_flag = 0 : say "Sorry it appears that $gc_file does not exist!";
		}
	}

	my ($n_file, $including_n);
	if(prompt_file("n"))
	{
		my $n_valid_flag = 1;
		$including_n = 1;
		while($n_valid_flag)
		{
			$n_file = prompt("Please enter the fullpath name for you N file: ");
			file_exists($n_file) ? $n_valid_flag = 0 : say "Sorry it appears that $n_file does not exist!";
		}
	}
	
	my ($map_file, $including_map);
	if(prompt_file("mappability"))
	{
		my $map_valid_flag = 1;
		$including_n = 1;
		while($map_valid_flag)
		{
			$map_file = prompt("Please enter the fullpath name for you mappability file: ");
			file_exists($map_file) ? $map_valid_flag = 0 : say "Sorry it appears that $map_file does not exist!";
		}
	}
}

########
# Subs #
########

# Prompt Yes or No for a file type
sub prompt_file
{
	my ($filetype) = @_;
	my $valid_flag = 1;
	my $answer;
	while($valid_flag)
	{
		$answer = uc(prompt("Include a $filetype file (Y/N)?  "));
		if($answer =~ m/^[YN]$/)
		{
			return($answer eq "Y");
		}
	}
}

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

# Operation Mode
sub handle_mode
{
	my ($name, $value) = @_;
	if($value eq "interactive")
	{
		$auto = 0;
	} elsif ($value eq "auto") {
		$auto = 1;
	} else {
		say "Error: Invalid value $value for mode arguement";
		say "Valid options are auto or interactive";
		die "$useage";
	}
}

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

# Generate the mosaicsFit command
# 1; analysis type
# 2: name of the bin
sub generate_fit_command 
{
	my ($analysis_type, $bin) = @_;
	my $fit = "fit".$analysis_type;
	my $fit_command = "$fit <- mosaicsFit($bin, analysisType=\"$analysis_type\")";
	return($fit_command);
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