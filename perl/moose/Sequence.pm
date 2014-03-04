package Sequence;
use Moose;
use namespace::autoclean;

has 'length' => (
	is  => 'ro',
	isa => 'Int',
	required => 1,
);

has 'seq' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

sub seq_array
{
	my $self = shift;
	my @seq = split(//, $self->seq);
	return \@seq;
}

__PACKAGE__->meta->make_immutable;
1;