#!/usr/bin/perl
# Parse ugly vcf into a nice table
use warnings;
use strict;
use feature qw(say);

my $dir = shift;
opendir DIR, $dir;
while(readdir(DIR))
{
	my $vcf = $_;
	next if($vcf !~ /\.vcf$/);
	my $out = $dir.$vcf.".out.txt";

	open VCF, "<", $dir.$vcf;
	open OUT, ">", $out;
	my @lines = grep {!/^#/} <VCF>;
	close VCF;
	say OUT $vcf;
	say OUT "LOC ID CONTROL EXPERIMENTAL";
	foreach my $line (@lines)
	{
		my @fields = split(/\t/,$line);
		my $loc = $fields[1];
		my $id = $fields[2];
		my $control = $fields[3];
		my $exp = $fields[4];
		say OUT "$loc $id $control --> $exp";
	}
	close OUT;
}
closedir(DIR);
