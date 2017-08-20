package Peak;
use Mouse;
use Number::Range;
use namespace::autoclean;

has ['peak_start','peak_stop'] => (
	isa => 'Int',
	required => 1,
	is => 'ro'
);

has 'peak_summit' => (
	isa => 'Int',
	is => 'ro',
	writer => '_set_summit'
);

has 'gtf' => (
	isa => 'Str',
	required => 1,
	is => 'ro'
);

has 'file' => (
	isa => 'Str',
	required => 1,
	is => 'ro'
);

has '_range' => (
	isa => 'Number::Range',
	is  => 'ro',
	writer => '_set_range'
);

sub BUILD {
	my $self = shift;

	if($self->peak_start > $self->peak_stop) {
		die "Cannot create peak with a start greater then the stop";
	}

	$self->_set_range(Number::Range->new($self->peak_start."..".$self->peak_stop));

	if($self->peak_summit ) {
		unless($self->_range->inrange($self->peak_summit)) {
			die "Cannot create peak with a summit not in range of start stop";
		}
	} else {
		my $peak_summit = int($self->peak_start+$self->peak_stop/2);
		$self->_set_summit($peak_summit);
	}
}

__PACKAGE__->meta->make_immutable;
1;