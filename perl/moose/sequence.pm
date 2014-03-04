package Sequence;
use Moose;
use namespace::autoclean

has 'length' => (
	is  => 'ro',
	isa => 'Int'
);

has 'seq' => (
	is => 'ro',
	isa => 'Str'
);



__PACKAGE__->meta->make_immutable;