use Genome;

my $fasta_file = shift;

my $genome  = Genome->new(
	fasta => $fasta_file,
);