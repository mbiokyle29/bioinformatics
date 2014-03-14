package Cigar;
use Moose;
use namespace::autoclean;

has 'raw_string' => (
	isa => 'Str',
	required => 1,
	is => 'ro'
);

has 'array' => (
	isa => 'ArrayRef',
	is => 'rw',
	auto_deref => 1,
	builder => '_build_array'
);

has 'stack' => (
	isa => 'Str',
	builder => '_build_stack',
	is => 'rw'
);

has 'length' => (
	isa => 'Int',
	builder => '_build_length',
	is => 'rw',
	lazy    => 1
);

has 'start_pos' => (
	isa => 'Int',
	builder => '_build_start',
	is => 'rw',
	lazy    => 1
);

has 'end_pos' => (
	isa => 'Int',
	builder => '_build_end',
	is => 'rw',
	lazy    => 1
);

sub _build_array
{
	my $self = shift;
	my $raw_string = $self->raw_string;
	my @array;

	while($raw_string =~ m/(\d+[SDMI])/g)
	{
		push(@array, $1);
	}

}

sub _build_stack
{
	my $self = shift;
	my $cigar = $self->raw_string;
	my @chunks;

	while($cigar =~ m/(\d+[SDMI])/g)
	{
		push(@chunks, $1);
	}

	my $cigar_stack;
	foreach my $chunk (@chunks)
	{
		my $code = chop($chunk);
		my $push = $code x $chunk;
		$cigar_stack.= $push;
	}
	return $cigar_stack;
}

sub _build_length
{
	my $self = shift;
	my $cigar = $self->raw_string;
	my @chunks;

	while($cigar =~ m/(\d+[DM])/g)
	{
		push(@chunks, $1);
	}

	my $cigar_length = 0;

	foreach my $chunk (@chunks)
	{
		if(chop($chunk) =~ m/[MD]/)
		{
			$chunk =~ m/(\d+)/;
			$cigar_length += $chunk;
		}
	}
	return $cigar_length;
}

sub _build_start
{
	return 1;
}

sub _build_end
{
	return 1;
}


1;
