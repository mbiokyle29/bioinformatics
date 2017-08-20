#
# Kyle McChesney
# File.pm
# 
package File;
use Mouse;
use File::Basename;
use File::Slurp;
use feature 'say';

use constant OK => 0;

has 'filename' => (
	isa => 'Str',
	is  => 'rw',
	required => 1
);

has 'fullpath' => (
	isa => 'Str',
	is  => 'rw'
);

has 'type' => (
	isa => 'Str',
	is  => 'rw',
	builder => '_build_type',
	lazy => 1
);

sub _build_type {
	my $self = shift;
	my $type = basename($self->fullpath);
	$type ? return $type : return -1; 
}

__PACKAGE__->meta->make_immutable();
1;