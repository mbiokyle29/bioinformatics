package Gene;
use Mouse;
use namespace::autoclean;

has ['gene_start','gene_stop'] => (
	isa => 'Int',
	required => 1,
	is => 'ro'
);

has 'name' => (isa => 'Str', required => 1, is => 'ro');
has 'gtf'  => (isa => 'Str', required => 1, is => 'ro');

sub BUILD {
	my $self = shift;
	if($self->gene_start < 0 || $self->gene_stop < 0) {
		die "Cannot create gene record with negative start stop";
	}
}

__PACKAGE__->meta->make_immutable;
1;