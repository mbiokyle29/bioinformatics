use Genome;

my $fasta_file = shift;
my $genome  = Genome->(
	fasta => $fasta_file
);