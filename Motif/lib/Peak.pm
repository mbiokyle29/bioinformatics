package Peak;
use Moose;
use Sequence;
use namespace::autoclean;
use Number::Range;
use Data::Printer;
use feature qw|say|;

has 'chromosome' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has id => (
	is => 'ro',
	isa => 'Str',
);

has ['start', 'summit', 'stop'] => (
	is => 'rw',
	isa => 'Int',
	required => 1,
);

has 'sequence' => (
	is => 'rw',
	isa => 'Object',
	required => 1,
);

has 'cobound' => (
	is => 'ro',
	isa => 'HashRef',
	required => 1,
);

# To handle the converstion of a string sequence to a seq object
around BUILDARGS => sub {
	my ($orig, $class, %args) = @_;
	
	my $seq_raw = $args{sequence};
	
	my $obj = Sequence->new (
		seq => $seq_raw
	);

	$args{sequence} = $obj;

	return $class->$orig(%args);
};

sub overlaps {
	my ($self, $left, $right) = @_;

	# Inverted bound check!
	unless(($self->stop < $left) or ($self->start > $right))
	{
		return 1;
	}

	return 0;
}

sub overlaps_lazy {
	my ($self, $other) = @_;
	my $left = $self->start - 100;
	my $right = $self->stop - 100;
	my $self_range = Number::Range->new($left."..".$right);

	if
	(
		$self_range->inrange($other->start)
		or $self_range->inrange($other->stop)
		or $self_range->inrange($other->summit)
	) 
	{ return 1; }

	else { return 0; }
}

# shallow equals
# if start, stop are the same
sub equals {
	my ($self, $other) = @_;

	unless($self->start == $other->start)
	{
		return 0;
	}

	unless($self->stop == $other->stop)
	{
		return 0;
	}

	return 1;
}

# string report
sub report {
	my $self = shift;

	my $string = $self->chromosome."\t".$self->start."\t".$self->stop;
	my $ref = $self->cobound;
	foreach my $prot (keys(%$ref))
	{
		if($$ref{$prot} == 1)
		{
			$string .= "\t".$prot;
		}
	}

	$string .= "\t".$self->id;
	$string .= "\n";

	return $string;
}

__PACKAGE__->meta->make_immutable;
1;