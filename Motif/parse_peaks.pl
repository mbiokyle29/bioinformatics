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

my ($username, $password, $do_insert);
GetOptions (
	"file=s" => \@peak_files,
	"uname=s" => \$username,
	"password=s" => \$password,
	"insert" => \$do_insert,
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
		&combine_peaks($matches, $dbh);
	}
}

sub combine_peaks {
	my $matches = shift;
	my $dbh = shift;

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
	my $end = $summit + 249;
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

	my $insert = "INSERT INTO refined_peaks (chromosome, start, stop, summit, sequence";
	my $values = " VALUES (?, ?, ?, ?, ?";
	foreach my $peak (@$matches)
	{
		my $ref_cobound = $peak->cobound;
		foreach my $protein (keys(%$ref_cobound))
		{
			$$ref_cobound{$protein} == 1 ? $cobound{$protein} = 1 : 0;
		}
	}

	my $cobound_count = 0;
	my @cobounds;
	for(grep {$cobound{$_} == 1} keys(%cobound))
	{
		push(@cobounds, 1);
		$insert .= ", ".$_;
		$values .= ", ?"; 
	}
	$insert = $insert.")".$values.")";
	my $sth = $dbh->prepare($insert);
	$sth->execute($chr, $start, $end, $summit, $substring, @cobounds);
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

sub clear_peak {
	my ($peak, $dbh) = @_;
	my $protein;
	
	# Cleared peak only has one cobound
	my $cobound = $peak->cobound;
	foreach my $key (keys(%$cobound))
	{
		if($$cobound{$key} == 1) { $protein .= $key; last; }
	}

	# Fix the length of the peak
	my $seq = $peak->sequence->seq;
	my $switch = 1;

	# Chomp first or last until get 499
	while(length($seq) > 499)
	{
		if($switch)
		{
			$seq = substr($seq, 1);
		} else {
			chop($seq);
		}
		$switch = $switch * (-1);
	}

	my $insert = "INSERT INTO refined_peaks (chromosome, start, stop, summit, sequence, $protein) VALUES (?, ?, ?, ?, ?, ?)";
	my $sth = $dbh->prepare($insert);cd 
	
	$sth->execute($peak->chromosome, $peak->start, $peak->stop, $peak->summit, $seq, 1);
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
			&clear_peak($current, $dbh);
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