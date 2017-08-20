#!/usr/bin/perl
use warnings;
use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;
use Getopt::Long;
use DBI;
use Peak;

# Global Peak Hash
my %peaks;

# Get all the peak files
my @peak_files;
my ($username, $password);
GetOptions (
	"file=s" => \@peak_files,
	"uname=s" => \$username,
	"password=s" => \$password
);

# DBI connect
my $dbh = DBI->connect('DBI:mysql:Peaks', $username, $password);

my $num_files = scalar(@peak_files);

# Iterate through all the peak files
foreach my $peak_file (@peak_files)
{
	my @lines = read_file($peak_file);

	# Bulid all the peak objects add to global hash
	foreach my $line (@lines)
	{
		next if ($line =~ /NA/);
		my $prot_name = $peak_file;
		$prot_name =~ s|^.*/||;
		my $peak = &build_peak($line, $prot_name, $dbh);
		
		if(exists($peaks{$peak->chromosome}))
		{

			push($peaks{$peak->chromosome}{peaks}, $peak);
			
			if($peak->start < $peaks{$peak->chromosome}{start})
			{
				$peaks{$peak->chromosome}{start} = $peak->start;
			}

			if($peak->stop > $peaks{$peak->chromosome}{stop})
			{
				$peaks{$peak->chromosome}{stop} = $peak->stop;
			}

		} else {
			$peaks{$peak->chromosome} = {
				start => $peak->start,
				stop => $peak->stop,
				peaks => [$peak]
			};
		}
	}
}

# We have all the peaks now
&bin_peaks_bystart(\%peaks, $num_files);


######################################################################
sub build_peak {
	my ($line, $peak_file, $keys) = @_;
	my @fields = split(/\t/, $line);

	# Get the scalar data members from array
	my ($chr, $start, $stop, $summit, $seq) = @fields[0,1,2,3,-1];

	# Build the dictionary
	my %dict = map { $_ => 0 } @$keys;
	$dict{$peak_file} = 1;


	# Create the peak object
	my $peak = Peak->new (
		chromosome => $chr,
		start => $start,
		summit => $summit,
		stop => $stop,
		sequence => $seq,
		prot => $peak_file,
		cobound => \%dict
	);
	p($peak);
	return $peak;
}

sub bin_peaks_bystart {
	my ($peaks, $num_files) = @_;
	# Each key of the peaks hash contains all the
	# info for each chromosome
	# so the peaks that fall in that chromosome and the min/max coordinates
	# that the peaks fall in  (for binning)
	foreach my $chr (keys(%$peaks))
	{
		my $min = $$peaks{$chr}{start};
		my $max = $$peaks{$chr}{stop};
		my $binsize = 10;

		# number of bins is range / binsize
		# extra points for hacked round up
		my $num_bins = $max - $min;
		$num_bins = $num_bins / $binsize;
		$num_bins = int(++$num_bins);

		say "min: $min   max: $max   bins: $num_bins for chr: $chr";
		my @bins;

		# some real computer science right here
		foreach my $peak (@{$$peaks{$chr}{peaks}})
		{
			# For an array of bins
			# the index an elment goes in
			# is equal to the elements value - starting value
			# 	integer divided by the binssize
			my $bindex = $peak->start - $min;
			$bindex = int($bindex / $binsize);
			if($bins[$bindex])
			{
				push($bins[$bindex], $peak);		
			} else {
				$bins[$bindex] = [$peak];
			}
		}

		# Take a look at the bins
		my $total = 0;
		my $matched = 0;
		my @unmatched;
		foreach my $bin (@bins)
		{
			next unless($bin);
			$total++;
			if(scalar(@$bin) == $num_files)
			{
				my $cobounds = $$bin[0]->cobound;
				foreach my $peak (@$bin)
				{
					if(Compare($peak->cobound, $cobounds))
					{
						$matched++;
					}
				}
			}
		}
		say "matched: $matched of $total";
	}
}

# Lets try string matching i guess?
