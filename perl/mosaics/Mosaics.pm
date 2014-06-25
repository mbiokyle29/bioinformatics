package Mosaics;
use Moose;
use namespace::autoclean;
use Statistics::R;

# Analysis type for MOSAICS fit
use constant OS => "OS";
use constant TS => "TS";
use constant IO => "IO";


has 'r_con' => (
	is => 'ro',
	isa => 'Object',
);

has 'r_log' => (
	is => 'rw',
	isa => 'Str'
);

has ['chip_file', 'input_file', 'analysis_type', 'file_format', 'out_loc', 'chip_bin', 'input_bin'] => (is => 'rw', isa => 'Str', lazy => 1);
has ['fragment_size', 'bin_size' ] => (is => 'rw', isa => 'Int', default => 200);
has ['map_score', 'gc_score', 'n_score'] => (is => 'rw', isa => 'Str', lazy => 1);
has ['bin_data'] => (is => 'rw', isa => 'Str', lazy => 1);

before 'make_chip_bin' => \&_have_file('chip', 'constructBins');
before 'make_input_bin' => \&_have_file('input', 'constructBins');
before 'make_chip_wiggle' => \&_have_file('chip', 'generateWig');
before 'make_input_wiggle' => \&_have_file('input', 'generateWig');
before 'read_bins' => \&_can_read_bins();

sub BUILD
{
	my $self = shift;
	$self->r_con = Statistics::R->new();
	$self->r_log = "R Connection Initialized";
	$self->_run_updates();
	$self->_load_libs();
}

### Sub for constructing a chip Bin
### Returns the name of the new chip bin (and overwrites the data memeber)
### 	or negative 1 if it fails
sub make_chip_bin
{
	my $self = shift; 
	my $const_bin = "constructBins(infile=\"".$self->chip_file."\", fileFormat=\"".$self->file_format."\", outfileLoc=\"".$self->out_loc."\", byChr=FALSE, fragLen=.".$self->fragment_length.", binSize=".$self->bin_size.")";
	if($self->r_con->run($const_bin))
	{
		$self->_log_command($const_bin);
		$chip_bin = $self->chip_file."_fragL".$self->fragment_length."_bin".$self->bin_size.".txt";
		$self->chip_bin($chip_bin);
		return $chip_bin;
	} else { return -1; }
}

### Sub for constructing a chip Bin
### Returns the name of the new input bin (and overwrites the data memeber)
### 	or negative 1 if it fails
sub make_input_bin
{
	my $self = shift; 
	my $const_bin = "constructBins(infile=\"".$self->input_file."\", fileFormat=\"".$self->file_format."\", outfileLoc=\"".$self->out_loc."\", byChr=FALSE, fragLen=.".$self->fragment_length.", binSize=".$self->bin_size.")";
	if($self->r_con->run($const_bin))
	{
		$self->_log_command($const_bin);
		$input_bin = $self->input_file."_fragL".$self->fragment_length."_bin".$self->bin_size.".txt";
		$self->input_bin($input_bin);
		return $input_bin;
	} else { return -1; }
}

### Generate chip wiggle file!
sub make_chip_wiggle
{
	my $self = shift;
	my $wiggle_command = "generateWig( infile=\"".$self->chip_file."\", fileFormat=\"".$self->in_format."\", outfileLoc=\"".$self->out_loc.")";
	$self->r_con->run($wiggle_command);
	$self->_log_command($wiggle_command);
}

### Generate input wiggle file!
sub make_input_wiggle
{
	my $self = shift;
	my $wiggle_command = "generateWig( infile=\"".$self->input_file."\", fileFormat=\"".$self->in_format."\", outfileLoc=\"".$self->out_loc.")";
	$self->r_con->run($wiggle_command);
	$self->_log_command($wiggle_command);
}

### Read in Bins
### WARNING AUTOMATICALLY USING SET DATA FIELDS! 
sub read_bins
{
	my $self = shift;
	$self->chip_file =~ m/^([\w]+)\..*$/;
	my $bin_data_name = $1.$self->analysis_type;
	
	# Set up strings for appending chosen data
	my $read_command = "$bin_data_name <- readBins";
	my $type_string = $self->_readbin_type_string();
	my $files_string = $self->_readbin_file_string();
	$read_command .= $type_string.", ".$files_string.")";
	
	$self->bin_data($bin_data_name);

}

## Compare set data memebers to analysis type
## Varify we can safely run readBins in MOSAICS
sub _can_read_bins
{
	my $self = shift;
	unless($self->analysis_type) { die "Cannot read bins without analysis_type being set!"; }
	unless($self->chip_bin)      { die "Cannot read bins without chip_bin file being set!"; }
	given($self->analysis_type)
	{
		when(OS) 
		{
			# Needs M + GC + N for OS 
			unless($self->map_score and $self->gc_score and $self->n_score) { 
				die "Cannot read bins in OS (one sample) mode without GC+M+N score incorporated!";
			}
		}
		when(TS or IO) {
			unless($self->input_bin and ) { die "Cannot read bins in two sample mode without input bin set"; }
		}
	}
}

# Internal validation Methods
## Sub for constructing an input Bin
sub _have_file
{
	my $self = shift;
	my $type = shift;
	my $operation = shift;

	unless($self->file_format) { die "Cannot perform $operation without a file format set!"; }
	if($type eq 'chip')
	{
		unless($self->chip_file) {
			die "Cannot perform $operation without a chip_file, please initialize!";
		}
	}
	elsif($type eq 'input')
	{
		unless($self->input_file) {
			die "Cannot perform $operation an input_file, please initialize!";
		}
	}	 
}

## R Library functions ##
sub _run_updates
{
	my $self = shift;
	my $connect = 'source("http://bioconductor.org/biocLite.R")';
	my $upgrader = 'biocLite()';
	my @commands = ($connect, $upgrader);
	$self->r_con->run(@commands);
	$self->_log_command($_) for @commands;
}

sub _load_libs
{
	my $self = shift;
	my $parallel = "library(parallel)";
	my $mosaics = "library(mosaics)";
	my @commands = ($parallel, $mosaics);
	$self->r_con->run(@commands);
	$self->_log_command($_) for @commands;
}

### Logging function to track all R commands ran
sub _log_command
{
	my ($self, $entry) = @_;
	$self->r_log .= $entry;
}

sub _readbin_type_string
{

}

sub _readbin_file_string
{

}