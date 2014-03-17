package Sam;
use Moose;
use Cigar;
use Sequence;
use namespace::autoclean;

has 'raw_string' => (
	isa => 'String',
	required => 1,
	is => 'ro'
);

has 'sequence' => (
	isa => 'Sequence',
	is => 'rw',
	builder => 
);

has 'cigar' => (
	isa => 'Cigar',
	is => 'rw'
);

sub BUILD
{
	my $self = shift;
}


__PACKAGE__->meta->make_immutable;
1;