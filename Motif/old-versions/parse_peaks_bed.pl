#!/usr/bin/perl
use warnings; use strict;
use feature qw|say|;
use File::Slurp;
use Data::Printer;
use Getopt::Long;
use DBI;
use lib "lib/"; use Peak;
use lib "/home/kyle/lab/Binder/lib"; use Binder;
use Data::Dumper;
use List::Util qw|sum min max|;


# Get all the peak files
my @peak_files;

our $cleared_bed = "cleared.bed";
our $merged_bed = "merged.bed";
our $misfits_bed = "misfits.bed";
our $combined_bed ="combined.bed";

write_file($cleared_bed, "track name=cleared description=\"cleared\"\n");
write_file($merged_bed, "track name=merged description=\"merged\"\n");
write_file($misfits_bed, "track name=misfits description=\"misfits\"\n");
write_file($combined_bed, "track name=combined description=\"Combined peaks\"\n");

my ($username, $password, $do_insert, $window_size);
GetOptions (
	"file=s" => \@peak_files,
	"uname=s" => \$username,
	"password=s" => \$password,
	"insert" => \$do_insert,
	"window=s" => \$window_size,
);

# DBI connect
my $dbh = DBI->connect('DBI:mysql:Peaks', $username, $password);

# FA dir

# Optional insert -- testing tsk tsk
if($do_insert) { &do_insert(\@peak_files, $dbh); }

# We have all the peaks in mysql now
my $chromosomes = &get_chromosomes($dbh);

# The magic
foreach my $chromosome (@$chromosomes){
	&iterative_search($dbh, $chromosome);
}

my @files = ($cleared_bed, $misfits_bed, $merged_bed);
&print_html(\@files);

######################################################################

#
# Peak Refinement - The magic method
#
sub refine_peak {
	my ($dbh, $matches) = @_;

	# matches is a group of overlapping peaks
	# determine if there is a shared point among
	# all the peaks
	
	# Sort by length low to high
	my @peaks_sorted = @$matches;
	my $all_overlap = &merge_overlap(\@peaks_sorted);
	
	# If yes then merge peak
	if($all_overlap)
	{
		append_file($merged_bed, $_->report) for @$matches;
		append_file($merged_bed, "\n");
		&combine_peaks($matches);
	}

	# If no then dump to misfits
	else {
		append_file($misfits_bed, $_->report) for @$matches;
		append_file($misfits_bed, "\n");
	}
}

sub combine_peaks {
	my $matches = shift;

	# Get weighted average summit point
	my $summit_sum = 0;
	$summit_sum += $_->summit for @$matches;
	
	my $denom = scalar(@$matches);

	my $summit = int(++$summit_sum/$denom);
	my $chr = $$matches[0]->chromosome;
	my $start = $summit - 250;

	# Peak at the left end of the chromo case
	if($start < 1) { $start = 1; }

	my $start_index = $start - 1;
	my $end = $summit + 250;
	my $length = $end - $start;

	# Get the sequence file
	my $fa_dir = "/data/k/web/";
	my $file = $fa_dir.$chr.".fa";
	my @sequence = read_file($file) or die "Cannot read";
	
	# clear the header line
	shift(@sequence);
	my $line = shift(@sequence);
	my $substring = substr($line, $start_index, $length);

	my %cobound = 
	(
		ebna2   => 0,
		ebna3a  => 0,
		ebna3b  => 0,
		ebna3c  => 0,
		rbpj92  => 0,
		rbpj234 => 0
	);

	foreach my $peak (@$matches)
	{
		my $ref_cobound = $peak->cobound;
		foreach my $protein (keys(%$ref_cobound))
		{
			$$ref_cobound{$protein} == 1 ? $cobound{$protein} = 1 : 0;
		}
	}

	my $peak_entry = "$chr\t$start\t$end\t";

	$peak_entry .= "-".$_ for(grep {$cobound{$_} == 1} keys(%cobound));
	append_file($combined_bed, $peak_entry."\n\n");

	say ">$chr $start $end\n$substring\n";
}


sub merge_overlap
{
	my $peak_stack = shift;
	
	# Base Case - last one in the stack
	if(scalar(@$peak_stack) == 1)
	{
		return 1;
	}

	# If not grab the first and second item pointers
	my $first = $$peak_stack[0];
	my $second = $$peak_stack[1];

	# Base Case - Test if they overlap
	unless($first->overlaps($second->start, $second->stop))
	{
		return 0;
	}

	# Merge and recurse
	else {	
		my $start_point = max($first->start, $second->start);
		my $end_point = min($first->stop, $second->stop);
		
		shift(@$peak_stack);
		
		# Create a new temporary peak which 
		# has a range that consists of the merged overlap range of
		# the first and second
		my $new = Peak->new (
			id => $$peak_stack[0]->id,
			start => $start_point,
			stop => $end_point,
			summit => $$peak_stack[0]->summit,
			chromosome => $$peak_stack[0]->chromosome,
			cobound => $$peak_stack[0]->cobound,
			sequence => $$peak_stack[0]->sequence->seq
		);

		# Overwrite the first peak in the stack
		$$peak_stack[0] = $new;

		# Recurse
		return &merge_overlap($peak_stack);
	}
}

sub clear_peak
{
	my $peak = shift;
	append_file($cleared_bed, $peak->report."\n");
}

#
# build peak - mini orm method
#
sub build_peak {
	my ($row_ref) = @_;
	
	my %cobound = (
		ebna2   => $$row_ref{ebna2  },
		ebna3a  => $$row_ref{ebna3a },
		ebna3b  => $$row_ref{ebna3b },
		ebna3c  => $$row_ref{ebna3c },
		rbpj92  => $$row_ref{rbpj92 },
		rbpj234 => $$row_ref{rbpj234},
	);

	my $peak = Peak->new (
		id         => $$row_ref{id},
		chromosome => $$row_ref{chromosome},
		start      => $$row_ref{start},
		stop       => $$row_ref{stop},
		summit     => $$row_ref{summit},
		sequence   => $$row_ref{sequence},
		cobound    => \%cobound
	);

	return $peak;	
}

#
# Fetching methods
#
sub fetch_peak {
	my ($dbh, $id) = @_;
	my $select = "SELECT * FROM peaks WHERE id=?";
	my $sth = $dbh->prepare($select);
	$sth->execute($id);

	my $row = $sth->fetchrow_hashref();
	return &build_peak($row);
}

sub fetch_all_peaks {
	my ($dbh, $chr) = @_;
	my $select = "SELECT * FROM peaks WHERE chromosome=? ORDER BY start ASC";
	my $sth = $dbh->prepare($select);
	$sth->execute($chr);

	my @peaks;

	while( my $row = $sth->fetchrow_hashref())
	{
		my $peak = &build_peak($row);
		push(@peaks, $peak);
	}

	return \@peaks;
}


#
# Inserting methods 
#
sub do_insert {
	my ($peak_files, $dbh) = @_;
	my $num_files = scalar(@$peak_files);

	# Iterate through all the peak files
	foreach my $peak_file (@$peak_files)
	{
		my @lines = read_file($peak_file);
	
		# Bulid all the peak objects add to global hash
		foreach my $line (@lines)
		{
			next if ($line =~ /NA/);
			my $prot_name = $peak_file;
			$prot_name =~ s|^.*/||;
			my $peak = &insert_peak($line, $prot_name, $dbh);
		}
	}
}

sub insert_peak {
	my ($line, $prot_name, $dbh) = @_;
	my @fields = split(/\t/, $line);

	# Get the scalar data members from array
	my ($chr, $start, $stop, $summit, $seq) = @fields[0,1,2,3,-1];

	# Set up query
	my $insert = "INSERT INTO peaks (chromosome, start, stop, summit, sequence, $prot_name) VALUES (?, ?, ?, ?, ?, ?);";
	my $one = 1;

	my $sth = $dbh->prepare($insert);
	$sth->execute($chr, $start, $stop, $summit, $seq, $one);

}


#
# Get various non peak stuff 
#
sub get_chromosomes {
	my $dbh = shift;
	my @chromosomes;
	
	my $query = "SELECT DISTINCT chromosome FROM peaks";
	my $sth = $dbh->prepare($query);
	$sth->execute;
	while(my $row = $sth->fetchrow_arrayref)
	{
		push(@chromosomes, $$row[0]);
	}
	return \@chromosomes;
}



sub iterative_search {
	my ($dbh, $chromosome) = @_;
	
	my $peaks = &fetch_all_peaks($dbh, $chromosome);

	# Find consecutive overlapping peaks
	# group them, and pass to refine when 
	# overlapping stops 
	my @group;
	my $current = shift(@$peaks);
	my $previous = $current;
	
	# The inital overlap boundry
	# For edge cases with long peaks starting 
	# before short, and then short not overlapping	
	my $group_left_boundry = $current->start;
	my $group_right_boundry = $current->stop;

	push(@group, $current);


	while(scalar(@$peaks))
	{
		# if a group falls in the overlap range
		if($$peaks[0]->overlaps($group_left_boundry, $group_right_boundry))
		{
			# Set it to current and push it to group
			$current = shift(@$peaks);
			push(@group, $current);

			# Update the group bounds
			if($current->start < $group_left_boundry)
			{
				$group_left_boundry = $current->start;
			}

			if($current->stop > $group_right_boundry)
			{
				$group_right_boundry = $current->stop;
			}

		# Check if this looks like a really alone peak
		} elsif (scalar(@group) == 1) {
			# handle the peak set that is finished
			&clear_peak($current);
			$previous = $current;
			
			# Set up the next current and the group
			$current = shift(@$peaks);
			$group_left_boundry = $current->start;
			$group_right_boundry = $current->stop;
			undef(@group);
			push(@group, $current);

		} else {
			# handle the peak set that is finished
			&refine_peak($dbh, \@group);
			$previous = $current;
			
			# Set up the next current and the group
			$current = shift(@$peaks);
			$group_left_boundry = $current->start;
			$group_right_boundry = $current->stop;
			undef(@group);
			push(@group, $current);
		}
	}
}


#
# Printing HTML files
#
sub print_html {
	my $files = shift;
	my $base = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=(?)%3A(?)-(?)";
	my $binder = Binder->new( base_string => $base);
	my $output = "peak_report.html";
	
	# Set up html header
	my $head = "<!DOCTYPE html>\n<html lang=\"en-US\">\n";
	$head .="<head>\n";
	$head .="\t<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap.min.css\" />\n";
	$head .="</head>\n<body>\n";
	write_file($output, $head);
	
	foreach my $file (@$files)
	{
		# set up the div
		my $div = "\t<div class=\"col-md-4\">\n\t\t<h1>$file</h1>\n";
		append_file($output, $div);

		# Read in the file and clear the track header
		my @lines = read_file($file);
		shift(@lines);
		
		my @starts;
		my @stops;
		my $chr;
		my @links;

		while(my $line = shift(@lines))
		{
			if($line eq "\n")
			{
				my $min = min(@starts);
				my $max = max(@stops);

				$min -= 100;
				$max += 100;

				$binder->bind($chr, $min, $max);
				push(@links, $binder->bound_string);
				undef(@starts);
				undef(@stops);
			}

			else
			{
				my ($t_chr, $t_start, $t_stop) = (split("\t", $line))[0..2];
				push(@starts, $t_start);
				push(@stops, $t_stop);
				$chr = $t_chr;
			}
		}
		my $out;
		my $counter = 1;
		foreach my $line (@links)
		{
			$out .= "\t\t<a href=\"$line\">$file - peak #$counter</a><br>\n";
			$counter++;
		}
		$out .= "\t</div>\n";
		append_file($output, $out);
	}
	append_file($output, "</body>\n");
}