#!/usr/bin/perl
use warnings;
use strict;
use DBI;
use Getopt::Long;
use File::Slurp;
use feature qw(say);


my ($genome, $username, $password);

GetOptions (
	"g=s" => \$genome,
	"u=s" => \$username,
	"p=s" => \$password,
);

unless ($genome && $username && $password)
{
	say "invalid command line args!"
	say "gen_genome.pl -g GENOME.fasta -u mysql-username -p mysql-password";
	die;
}
$genome =~ m/(.+)\.fa.+$/;
my $genome_name = $1;

my $dbh = DBI->connect('dbi:mysql:genomes',$username,$password);

# TODO check if that genome already exists
	# if so die

# TODO make table
my $make_table = "CREATE TABLE ? (base VARCHAR(1), position long)";

#TODO foreach base insert!!? bad bad oh well