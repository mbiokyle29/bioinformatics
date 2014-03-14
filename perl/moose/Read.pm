package Read;
use Moose;
use Cigar;
use Sequence;
use namespace::autoclean;

has 'raw_read' => (
	isa => 'String',
	required => 1,
	is => 'ro'
);

has 'sequence' => (
	isa => 'Object'
);

has 'cigar' => (
	isa => 'Object'
);

sub BUILD
{
	my $self = shift;
}


__PACKAGE__->meta->make_immutable;
1;