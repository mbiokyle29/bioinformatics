package Genome;
use Moose;
use namespace::autoclean;
use DBI;

my $dbh = DBI->connect('dbi:mysql:genomes','genomes','ebvHACK958$');

has 'length' => (
	is  => 'ro',
	isa => 'Int',
);

has 'seq' => (
	is => 'ro',
	isa => 'Str',
);

has 'fasta' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'name' => (
	is => 'ro',
	isa => 'Str',
);


sub build
{
	my $self = shift;
	my $fasta = $self->fasta;
	open my $fa, "<", $fasta;
	
	my $ref_name = <$fa>;
	$ref_name =~ s/\.fa(sta)*//;

	#Check if already exists save us some time
	my $table_check = "SHOW TABLES LIKE '?'";
	my $check_sth = $dbh->prepare($table_check);
	$check_sth->execute($ref_name);

	if($check_sth->fetchrow_arrayref->[0])
	{
		#no need to dump to mysql
		# jsut build model of existing table
		$self->name = $ref_name;
	}

	my $ref_seq;
	while(<$fa>)
	{
		my $line = $_;
		if($line =~ m/^[^>\#]/)
		{
			$ref_seq.=$line;
		}
	}
	my $make_table = "CREATE TABLE ? (base VARCHAR(1), position long)";
	my $sth = $dbh->prepare($make_table);
	$sth->execute($ref_name) or die "Error Creating table for genome $ref_name $DBI::errstr\n";	

	my @bases = split(//,$ref_seq);
	my $counter = 1;
	my $insert = "INSERT INTO $ref_name (base, position) VALUES(?, ?)";
	foreach my $base (@bases)
	{
		my $insert_h = $dbh->prepare($insert);
		$insert_h->execute($base, $counter);
		$counter++;
	}
}

__PACKAGE__->meta->make_immutable;
1;