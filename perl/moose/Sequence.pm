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

sub base_at
{
	my $self = shift;
	my $pos = shift;
	$pos--
;	return substr($self->seq, $pos,1); 
}

sub range
{
	my $self = shift;
	my $start = shift;
	$start--;
	my $end = shift;
	return substr($self->seq, $start, $end);
}

__PACKAGE__->meta->make_immutable;
1;