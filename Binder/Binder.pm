package Binder;
use Mouse;
use namespace::autoclean;
our $VERSION = "1.0";

has 'base_string' => (is => 'rw', isa => 'Str', required => 1);
has 'slot_number' => (is => 'rw', isa => 'Int');
has 'bound_string' => (is => 'rw', isa => 'Str');
has 'bind_symbol' => (is => 'rw', isa => 'Str', default => '(?)');

sub BUILD {
	my $self = shift;
	$self->set_slot_number();
}

sub bind {
	my $self = shift;
	my @params = @_;
	$self->set_slot_number();
	if(scalar(@params) != $self->slot_number)
	{
		die "Cannot bind values into base string, slot number doesn't match params given";
	}
	my $temp = $self->base_string;
	for my $param (@params)
	{
		my $bind = quotemeta($self->bind_symbol);
		$temp =~ s/$bind/$param/;
	}
	$self->bound_string($temp);
	return $self->bound_string();
}

sub set_slot_number 
{
	my $self = shift;
	my $bind_symbol = quotemeta($self->bind_symbol);
	my $slot_number = () = ($self->base_string =~ m/$bind_symbol/g);
	$self->slot_number($slot_number);
}

__PACKAGE__->meta->make_immutable;
1;
