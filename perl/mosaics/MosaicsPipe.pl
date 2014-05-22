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
	"format=s" => \&handle_chip,
	"input=s" => \$input,
	"map=s" => \$map_score,
	"gc=s" => \$gc_score,
	"n=s" => \$n_score
);
unless($analysis_type and $auto)
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
	my $bin_name = "bin1";
	my $bin_command = $bin_name." <- readBins(";
	my $chip_bin;
	my $input_bin;
	
	say "Starting automatic run-mode";
	say "Analysis Type => $analysis_type";
	
	# Handle the chip file
	$chip = abs_path($chip);
	my $out_path = $chip;
	$out_path =~ s/\/[^\/]+$/\//;
	my $const_chip = "constructBins(infile=\"$chip\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
	$chip_bin = $chip."_fragL200_bin50.txt";
	
	$type_string.='"chip", ';
	$files_string.="\"$chip_bin\", ";

	push(@commands, $const_chip);

	if($input)
	{
		$input = abs_path($input);
		my $const_input = "constructBins(infile=\"$input\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
		push(@commands, $const_input);
		$input_bin = $input."_fragL200_bin50.txt";
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

	chomp($type_string);
	chomp($files_string);
	
	$type_string =~ s/,\s$/)/;
	$files_string =~ s/,\s$/)/;
	$bin_command .= $type_string.", ".$files_string.")";
	push(@commands, $bin_command);
	
	# Set up R env
	say &run_updates($r_con);
	say &load_libs($r_con);
	say "\nLibraries loaded correctly being processing results \n";

	foreach my $command (@commands)
	{
		say "running $command";
		say $r_con->run($command);
	}
	say $r_con->run("show($bin_name)");

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
sub handle_chip 
{
	my ($name, $value) = @_;
	my @valid_types = qw|eland_result eland_extended eland_export bowtie sam bed csem|;
	unless($value ~~ @valid_types)
	{
		say "CHiP File format $value is invald!";
		say "Must be one of: @valid_types";
		die "$useage";
	}
	$chip_type = $value;
}