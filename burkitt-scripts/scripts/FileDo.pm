#
# Kyle McChesney
# FileDo.pm
# 
package FileDo;
use Mouse;
use File;
use File::Slurp;
use Data::Printer;
use File::Basename;
use Binder;
use feature 'say';

use constant OK => 1;

our $VERSION = "1.0";

# Array of file objects
has 'files' => (
	isa => 'ArrayRef[File]',
	is  => 'rw'
);

# Only required, dir to read things from
has 'basedir' => (
	isa => 'Str',
	is  => 'ro',
	required => 1
);

has 'filetype' => (
	isa => 'Str',
	is  => 'ro',
);

sub BUILD {
	my $self = shift;
	my @files = read_dir($self->basedir, prefix => 1);
	my @regex_matched;

	if($self->filetype) {
		my $filetype = $self->filetype;
		foreach my $file (@files) {
			if ($file =~ m/.*\.$filetype$/) {
				push(@regex_matched, $file);
			}
		}
		@files = @regex_matched;
	}

	my @file_objs;
	foreach my $file (@files) {
		my $file_name = fileparse($file);
		my $file = File->new(filename => $file_name, fullpath => $file);
		push(@file_objs, $file);
	}
	$self->files(\@file_objs);
}

sub print {
	my $self = shift;
	say "Dumping File list";
	foreach my $file (@{$self->files}) {
		say $file->filename();
	}
	return OK;
}

sub do_name {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);
	foreach my $file (@{$self->files}) {
		$binder->bind($file->filename);
		say " Running " .$binder->bound_string;
		system($binder->bound_string) == 0 or die $binder->bound_string." command failed!";
	}
	return OK;
}

sub print_name {
	my ($self, $code) = @_;
	my $return;
	my $binder = Binder->new(base_string => $code);
	foreach my $file (@{$self->files}) {
		$binder->bind($file->filename);
		$return .= $binder->bound_string();
	}
	return $return;
}

sub do_fp {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);
	foreach my $file (@{$self->files}) {
		$binder->bind($file->fullpath);
		say " Running " .$binder->bound_string;
		system($binder->bound_string) == 0 or die $binder->bound_string." command failed!";
	}
	return OK;
}

sub do_name_fp {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);
	foreach my $file (@{$self->files}) {
		$binder->bind($file->filename, $file->fullpath);
		say " Running " .$binder->bound_string;
		system($binder->bound_string) == 0 or die $binder->bound_string." command failed!";
	}
	return OK;
}

sub do_all {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);
	my $list;
	
	foreach my $file (@{$self->files}) {
		$list .= $file->fullpath.",";
	}
	$list =~ s/,$//;
	
	$binder->bind($list);
	say "RUNNING:";
	say $binder->bound_string;
	system($binder->bound_string) == 0 or return 0;
	return OK;
}

sub do_name_name {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);
	foreach my $file (@{$self->files}) {
		$binder->bind($file->filename, $file->filename);
		say " Running " .$binder->bound_string;
		system($binder->bound_string) == 0 or die $binder->bound_string." command failed!";
	}
	return OK;
}

sub do_mated_name_fp_fp {
	my ($self,$code) = @_;
	my $binder = Binder->new(base_string => $code);

	unless($binder->slot_number() == 3) {
		say " mated name fp fp requires a three slot input code";
		return -1;
	}

	my %pairs;
	foreach my $file (@{$self->files}) {
		if($file->filename =~ m/(.*)_\d\.fastq/) {
			push(@{ $pairs{$1} }, $file);
		} else {
			say "A file in the list was not in the correct mated form!";
			return -1;
		}
	}

	foreach my $file_base (keys(%pairs)) {
		p(@{$pairs{$file_base}});
		my ($file_one, $file_two) = sort {$a->filename cmp $b->filename} @{$pairs{$file_base}};
		$binder->bind($file_base, $file_one->fullpath, $file_two->fullpath );
		system($binder->bound_string) == 0 or die $binder->bound_string." failed";
	}

	return OK;
}

sub do_mated_basename_basename {
	my ($self,$code) = @_;
	my $binder = Binder->new(base_string => $code);

	unless($binder->slot_number() == 2) {
		say " mated basename basename requires a two slot input code";
		return -1;
	}

	my %pairs;
	foreach my $file (@{$self->files}) {
		if($file->filename =~ m/(.*)_\d\.fastq/) {
			push(@{ $pairs{$1} }, $file);
		} else {
			say "A file in the list was not in the correct mated form!";
			return -1;
		}
	}

	foreach my $file_base (keys(%pairs)) {
		$binder->bind($file_base,$file_base);
		system($binder->bound_string) == 0 or die $binder->bound_string." failed";
	}
	return OK;
}

sub do_basename {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);

	my %pairs;
	foreach my $file (@{$self->files}) {
		if($file->filename =~ m/(.*)_\d\.fastq/) {
			push(@{ $pairs{$1} }, $file);
		} else {
			say "A file in the list was not in the correct mated form!";
			return -1;
		}
	}

	foreach my $file_base (keys(%pairs)) {
		$binder->bind($file_base);
		system($binder->bound_string) == 0 or die $binder->bound_string." failed";
	}
	return OK;
}

sub print_basename {
	my ($self, $code) = @_;
	my $binder = Binder->new(base_string => $code);
	my $return;
	my %pairs;
	foreach my $file (@{$self->files}) {
		if($file->filename =~ m/(.*)_\d\.fastq/) {
			push(@{ $pairs{$1} }, $file);
		} else {
			say "A file in the list was not in the correct mated form!";
			return -1;
		}
	}

	foreach my $file_base (keys(%pairs)) {
		$binder->bind($file_base);
		$return .= $binder->bound_string == 0 or die $binder->bound_string." failed";
	}
	return $return;
}

__PACKAGE__->meta->make_immutable();
1;
