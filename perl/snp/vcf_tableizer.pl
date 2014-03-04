#!/usr/bin/perl
# Parse ugly vcf into a nice table
use warnings;
use strict;
use feature qw(say switch);
use Data::Dumper;

my $dir = shift;
opendir DIR, $dir;
while(readdir(DIR))
{
	my $vcf = $_;
	next if($vcf !~ /\.vcf$/);
	my $out = $dir.$vcf.".out.txt";
	
	open VCF, "<", $dir.$vcf;
	open OUT, ">", $out;
	
	my @header;
	my @vars;
	my @top;
	
	while(<VCF>)
	{
		chomp();
		my $line = $_;
		
		for($line) 
		{
			when (m/^##INFO/) {push(@header,$line)}
			when (m/^#\w/) {@top = split("\t",$line);} 
			when (m/^[^#]/) {push(@vars,$line)}
		}
	}
	
	my %infos;	
	foreach my $line (@header)
	{
		$line =~ m/INFO=<ID=(\w*),.+Description="(.*)">$/;
		$infos{$1} = $2;
	}
	
	foreach my $tops(@top)
	{
		if ($tops =~ m/INFO/) 
		{
			print OUT join("\t", sort(keys(%infos)));
			print OUT "\t";
		}
		else { print OUT $tops."\t"; }
	}
	say OUT "";
	
	foreach my $line (@vars)
	{
		
		my %values;
		foreach my $key (keys(%infos)) {$values{$key}="X";}
		my @split = split("\t",$line);
		for(my $i = 0; $i < @split; $i++)
		{
			if($i == 7)
			{
				my @info_values = split(";",$split[$i]);
				foreach my $info_pair (@info_values)
				{
					$info_pair =~ m/(.+)=(.+)/;
					$values{$1} = $2;
				}
				foreach my $key (sort(keys(%values)))
				{
					print OUT $values{$key}."\t";
				}
			}
			else { print OUT "$split[$i] \t"; }
		}
		say OUT "";
	}
	
}
closedir(DIR);
