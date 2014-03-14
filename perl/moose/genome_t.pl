use lib "/home/kyle/lab/perlpipe/perl/moose/";
use Genome;


my $fasta_file = shift;

my $genome  = Genome->new(
	fasta => $fasta_file,
);