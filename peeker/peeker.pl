#!/usr/bin/env perl
#
# Kyle McChesney
# peeker.pl
#  top level script for peeker web app
use Mojolicious::Lite;
use DBI;
use Data::Printer;
use lib "./lib";
use Peak;
use Crud;
use Gene;

our $config = plugin 'JSONConfig';
our $crud = Crud->new(config => $config);

get '/' => sub {
	my $self = shift;
	my $peaks = $crud->get_peaks();
	$self->stash( peaks => $peaks);
	my $genes = $crud->get_genes();
	$self->stash( genes => $genes);
	$self->render('index');
};

get '/peaks' => sub {
	my $self = shift;
	my $peaks = $crud->get_peaks();
	$self->stash( peaks => $peaks);
	$self->render();
};

get '/genes' => sub {
	my $self = shift;
	my $genes = $crud->get_genes();
	$self->stash( genes => $genes);
	$self->render('genes');
};

get '/query' => sub {
	my $self = shift;
	my $peaks = $crud->get_peaks();
	$self->stash( peaks => $peaks);
	my $genes = $crud->get_genes();
	$self->stash( genes => $genes);
	$self->render('query');
};

post '/upload_peaks' => sub {
	my $self = shift;
	my $file = $self->param('peak_file');
	&parse_peak_file($file);
	my $peaks = $crud->get_peaks();
	$self->stash( peaks => $peaks);
	$self->redirect_to('peaks');
};

post '/upload_genes' => sub {
	my $self = shift;
	my $file = $self->param('gene_file');
	&parse_gene_file($file);
	my $genes = $crud->get_genes();
	$self->stash( genes => $genes);
	$self->redirect_to('genes');
};


sub parse_gene_file {
	my $file = shift;
	my @lines = split("\n", $file->slurp());
	foreach my $peak (@lines) {
		my ($start, $stop, $name, $gtf) = split("\t",$peak);
		$crud->insert_gene($start, $stop, $name, $gtf);
	}
}

sub parse_peak_file {
	my $file = shift;
	my @lines = split("\n", $file->slurp());
	foreach my $peak (@lines) {
		my ($start, $stop, $summit, $gtf) = split("\t",$peak);
		$crud->insert_peak($start, $stop, $summit, $gtf, $file->filename);
	}
}

sub find_closest_above {
	my ($point, $canidates) = @_;
	
}

sub find_closest_below {
	
}

app->start;