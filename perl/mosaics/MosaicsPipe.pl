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
use feature qw|say|;
use File::Slurp;
use Statistics::R;

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

# Check we got everything for automatic mode
# Then run
if($auto && &verify_auto_reqs()) 
{ 
	my $type_string = "type = c(";
	my $files_string = "fileName = c(";
	my $bin_command = "bin1 <- readBins(";
	
	say "Starting automatic run-mode";
	say "Analysis Type => $type";
	
	# Handle the chip file
	my $out_path = abs_path($chip);
	my $const_chip = "constructBins(infile=\"$chip\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
	push(@commands, $const_chip);

	if($input)
	{
		my $const_input = "constructBins(infile=\"$input\", fileFormat=\"$chip_type\", 
		outfileLoc=\"$out_path\", byChr=FALSE, fragLen=200, binSize=50)";
		push(@commands, $const_input);
	}


	

	bin1 <- readBins(type = c("chip", "input", "M", "GC", "N"), 
		fileName = c("chip.sam_fragL200_bin50.txt", "input.sam_fragL200_bin50.txt", 
			"hg19_ebv_M_fragL200_bin50.txt", "hg19_ebv_GC_fragL200_bin50.txt", "hg19_ebv_N_fragL200_bin50.txt"))

}


































# Verify that user inputted the nessecary data for the given analysis type
sub verify_auto_reqs
{
	# At the very least we need the chip data no matter the type
	unless($chip) 
	{ 
		die "chip data must be submitted for automatic OS analysis \n $usage"; 
	}
	given($type)
	{
		when(OS)
		{
			unless($input or ($map_score and $gc_score and $n_score))
			{
				die "input data or GC,M and N data must be submitted for automatic OS analysis \n $usage";
			}
		}
		when(TS)
		{
			unless($input)
			{
				die "input data must be submitted for automatic TS analysis \n $usage";
			}
		}
		when(IO)
		{
			unless ($input) 
			{
				die "input data must be submitted for automatic IO analysis \n $usage";
			}
		}
	}
	return 1;
}

## Libraries ##
sub run_updates
{
	my $r_con = shift;
	my $connect = 'source("http://bioconductor.org/biocLite.R")';
	my $installer = 'biocLite("BiocInstaller")';
	my $upgrader = 'biocLite("BiocUpgrade")';
	my $upgrade = 'update.packages()';
	my @commands = ($connect, $installer, $upgrader, $upgrade);
	$r_con->run(@commands) or die "Could not check for updates";
	return 1;
}

sub load_libs
{
	my $r_con = shift;
	my $parallel = "library(parallel)";
	my $mosaics = "library(mosaics)";
	$r_con->run(($parallel, $mosaics)) or die "Could not load R libraries";
	return 1;
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
		die "$usage";
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
		die "$usage";
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
		die "$usage";
	} 
}