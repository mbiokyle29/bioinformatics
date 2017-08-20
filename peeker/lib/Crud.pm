package Crud;
use Mouse;
use namespace::autoclean;

has 'config' => (isa => 'HashRef', is => 'ro', required => 1);

sub _connect
{
	my $self = shift;
	my $dbh = DBI->connect($self->config->{dsn}, $self->config->{user}, $self->config->{pass}) or die "Could not connect to DB $DBI::errstr";
	return $dbh;
}

sub insert_peak {
	my ($self, $peak_start, $peak_stop, $peak_summit, $gtf, $file) = @_;
	my $dbh = $self->_connect();
	
	unless($peak_start && $peak_stop && $peak_summit && $gtf && $file) {
		die "Cannot create peak without all values";
	}

	my $peak_insert = "INSERT INTO peaks (peak_start,peak_stop,peak_summit,gtf,file) VALUES (?,?,?,?,?)";

	my $peaks_handle = $dbh->prepare($peak_insert);
	$peaks_handle->execute($peak_start, $peak_stop, $peak_summit, $gtf, $file) or die "Could not insert peak:  $DBI::errstr";
}

sub get_peaks {
	my $self = shift;
	my $dbh = $self->_connect();
	my $peaks_get = "SELECT * FROM peaks";
	my $peaks_h = $dbh->prepare($peaks_get);
	$peaks_h->execute();
	my %peaks;
	while(my $hash = $peaks_h->fetchrow_hashref()) {
		my $peak = Peak->new(
			peak_start  => $$hash{'peak_start'},
			peak_stop   => $$hash{'peak_stop'},
			peak_summit => $$hash{'peak_summit'},
			gtf         => $$hash{'gtf'},
			file        => $$hash{'file'}
		);
		push( @{$peaks{$$hash{'file'}}}, $peak );
	}
	return \%peaks;
}

sub get_genes {
	my $self = shift;
	my $dbh = $self->_connect();
	my $genes_get = "SELECT * FROM genes";
	my $genes_h = $dbh->prepare($genes_get);
	$genes_h->execute();
	my @genes;
	while(my $hash = $genes_h->fetchrow_hashref()) {
		my $gene = Gene->new(
			gene_start  => $$hash{'gene_start'},
			gene_stop   => $$hash{'gene_stop'},
			name 		=> $$hash{'name'},
			gtf         => $$hash{'gtf'},
		);
		push( @genes, $gene );
	}
	return \@genes;
}

sub get_gene {
	my ($self, $gene_id) = @_;
	my $dbh = $self->_connect();
	my $genes_get = "SELECT * FROM genes WHERE gene_id = ?";
	my $genes_h = $dbh->prepare($genes_get);
	$genes_h->execute($gene_id);
	my $hash = $genes_h->fetchrow_hashref();
	my $gene = Gene->new(
		gene_start  => $$hash{'gene_start'},
		gene_stop   => $$hash{'gene_stop'},
		name 		=> $$hash{'name'},
		gtf         => $$hash{'gtf'},
	);
	die "Couldnt fetch gene: $gene_id" unless defined($gene);
	return $gene;
}

sub insert_gene {
	my ($self, $gene_start, $gene_stop, $name, $gtf) = @_;
	my $dbh = $self->_connect();
	
	unless($gene_start && $gene_stop && $name && $gtf) {
		die "Cannot create gene without all values";
	}

	my $gene_insert = "INSERT INTO genes (gene_start,gene_stop,name,gtf) VALUES (?,?,?,?)";

	my $gene_handle = $dbh->prepare($gene_insert);
	$gene_handle->execute($gene_start, $gene_stop, $name, $gtf) or die "Could not insert peak:  $DBI::errstr";
}

__PACKAGE__->meta->make_immutable;
1;