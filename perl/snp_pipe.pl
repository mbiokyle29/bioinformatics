#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my ($ref, $bam_dir, $gatk);

GetOptions
(
	"d=s" => \$bam_dir,
	"m=s" => \$ref,
	"g=s" => \$gatk 
);

# Grab all the files, ignore . and ..
opendir DIR, $bam_dir;
my @bam_files = grep { !/^\./ } readdir(DIR);
closedir DIR;

for my $bam_file (@bam_files)
{	
	next unless($bam_file =~ m/^(.*)\.bam$/);
	$bam_file =~ m/(.*)\.bam$/;
	my $out = $1.".vcf";
	my $output = `java -d64 -Xms3g -Xmx4g -jar $gatk -T UnifiedGenotyper -R $ref -I $bam_file -o $out -nct 20 -glm BOTH -minIndelFrac 0.5 -baq CALCULATE_AS_NECESSARY -dcov 1000 -A AlleleBalance -A Coverage -A MappingQualityZero -stand_emit_conf 10.0 -rf BadCigar -rf NotPrimaryAlignment`;
	print $output;
	print "\n\n";
}
