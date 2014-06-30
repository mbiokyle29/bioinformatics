#!/usr/bin/perl
use warnings;
use strict;
use Statistics::R;
use Data::Printer;

use feature qw|say switch|;

# R connection
#my $r_con = Statistics::R->new();
#print"->"; say $r_con->run("library(mosaics)");
#print"->"; say $r_con->run("library(mosaicsExample)");
#print"->"; say $r_con->run("data(exampleFit)");
#print"->"; &try_fit("IO")

my @bgEst = qw|matchLow rMOM|;
my @signalModes = qw|1S 2S|;
my @fdrs = (.1, .08, .05);
my @truncProb = (.999, .95, .80);
p(@fdrs);



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

sub try_peak
{
	my ($fit, $r_con) = @_;
	my @fdr_range = (0.05, 0.08, 0.1);
	my @thres_range = 11..20;
	my @minsizes = (55..60, 90..95, 195..200);
	my @max_gaps = (200, -1);
	for my $fdr (@fdr_range)
	{
		for my $thres (@thres_range)
		{
			for my $minsize (@minsizes)
			{
				for my $max (@max_gaps)
				{
					# Call peaks
					eval {  };
				}
			}
		}
	}
}