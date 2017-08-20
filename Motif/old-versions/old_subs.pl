sub get_min {
	my ($dbh,$chr) = @_;

	my $query = "SELECT start FROM peaks WHERE chromosome=? ORDER BY start ASC LIMIT 1";
	my $sth = $dbh->prepare($query);
	$sth->execute($chr);

	my $res = @{$sth->fetchrow_arrayref()}[0];
	return $res;
}

sub get_max {
	my ($dbh,$chr) = @_;

	my $query = "SELECT stop FROM peaks WHERE chromosome=? ORDER BY stop DESC LIMIT 1";
	my $sth = $dbh->prepare($query);
	$sth->execute($chr);

	my $res = @{$sth->fetchrow_arrayref()}[0];
	return $res;
}

sub get_lengths {
	my ($dbh,$chr) = @_;
	my @lengths;

	my $query = "SELECT start, stop FROM peaks WHERE chromosome=?";
	my $sth = $dbh->prepare($query);
	$sth->execute($chr);


	while(my $row = $sth->fetchrow_arrayref())
	{
		my ($start, $stop) = @$row;
		my $length = $stop-$start;
		push(@lengths, $length);
	}
	return \@lengths;
}

sub get_average_length {
	my ($dbh, $chr) = @_;
	my $lengths = &get_lengths($dbh, $chr);
	my $average = sum(@$lengths)/scalar(@$lengths);
	say "$chr --> $average";
	return $average;
}

#
# Identical site methods
#
sub find_identical_associations {
	my ($dbh, $chromosome, $score_hash) = @_;

	# check peaks with identical starts
	# first get a list of positions that more then one
	# peaks start at
	my $start_list = &get_identical_starts($dbh, $chromosome);
	my @start_ids;

	# Then fetch all the peaks at each position
	# and update their scores
	foreach my $start (@$start_list)
	{
		my $ids = &fetch_identical_starts($dbh, $chromosome, $start);

		my @ids_without;

		foreach my $id (@$ids)
		{
			@ids_without = grep { $_ != $id } @$ids;
			$$score_hash{$id}{$_}{start_identical} = 1 for @ids_without;
		}
	}

	# Now do identical ends
	my $ends_list = &get_identical_ends($dbh, $chromosome);
	my @ends_ids;

	# Then fetch all the peaks at each position
	# and update their scores
	foreach my $end (@$ends_list)
	{
		my $ids = &fetch_identical_ends($dbh, $chromosome, $end);
		my @ids_without;

		foreach my $id (@$ids)
		{
			@ids_without = grep { $_ != $id } @$ids;
			$$score_hash{$id}{$_}{end_identical} = 1 for @ids_without;
		}
	}

	# Now do identical summits
	my $summits_list = &get_identical_summits($dbh, $chromosome);
	my @summits_ids;

	# Then fetch all the peaks at each position
	# and update their scores
	foreach my $summit (@$summits_list)
	{
		my $ids = &fetch_identical_summits($dbh, $chromosome, $summit);
		my @ids_without;

		foreach my $id (@$ids)
		{
			@ids_without = grep { $_ != $id } @$ids;
			$$score_hash{$id}{$_}{summit_identical} = 1 for @ids_without;
		}
	}
}

sub get_identical_starts {
	my ($dbh, $chr) = @_;
	my $list_query = "SELECT start FROM peaks WHERE chromosome = ? GROUP BY start HAVING COUNT(*) > 1";
	
	my $sth = $dbh->prepare($list_query);
	$sth->execute($chr);

	my @identicals;
	while(my $row = $sth->fetchrow_arrayref())
	{
		push(@identicals, $$row[0]);
	}
	return \@identicals;
}

sub fetch_identical_starts {
	my ($dbh, $chr, $iden) = @_;
	my $get_query = "SELECT id FROM peaks WHERE chromosome=? AND start=?";
	my $sth = $dbh->prepare($get_query);
	$sth->execute($chr, $iden);

	my @ids;
	while( my $row = $sth->fetchrow_arrayref())
	{
		push (@ids, $$row[0]);
	}

	return \@ids;
}

sub get_identical_ends {
	my ($dbh, $chr) = @_;
	my $list_query = "SELECT stop FROM peaks WHERE chromosome = ? GROUP BY stop HAVING COUNT(*) > 1";
	
	my $sth = $dbh->prepare($list_query);
	$sth->execute($chr);

	my @identicals;
	while(my $row = $sth->fetchrow_arrayref())
	{
		push(@identicals, $$row[0]);
	}
	return \@identicals;
}

sub fetch_identical_ends{
	my ($dbh, $chr, $iden) = @_;
	my $get_query = "SELECT id FROM peaks WHERE chromosome=? AND stop=?";
	my $sth = $dbh->prepare($get_query);
	$sth->execute($chr, $iden);

	my @ids;
	while( my $row = $sth->fetchrow_arrayref())
	{
		push (@ids, $$row[0]);
	}

	return \@ids;
}

sub get_identical_summits {
	my ($dbh, $chr) = @_;
	my $list_query = "SELECT summit FROM peaks WHERE chromosome = ? GROUP BY summit HAVING COUNT(*) > 1";
	
	my $sth = $dbh->prepare($list_query);
	$sth->execute($chr);

	my @identicals;
	while(my $row = $sth->fetchrow_arrayref())
	{
		push(@identicals, $$row[0]);
	}
	return \@identicals;
}

sub fetch_identical_summits { 
	my ($dbh, $chr, $iden) = @_;
	my $get_query = "SELECT id FROM peaks WHERE chromosome=? AND summit=?";
	my $sth = $dbh->prepare($get_query);
	$sth->execute($chr, $iden);

	my @ids;
	while( my $row = $sth->fetchrow_arrayref())
	{
		push (@ids, $$row[0]);
	}

	return \@ids;
}

#
# Close site methods
#
sub find_close_associations 
{
	my ($dbh, $chr, $min, $max, $window_size, $score_hash) = @_;

	# Set up inital window
	my $left = $min;
	my $right = $min+$window_size;
	my $increment = int(++$window_size/2);

	# main loop
	while($right <= $max)
	{
		# For holding temp list
		my @ids_without;

		my $starts = &starts_in_window($dbh, $chr, $left, $right);
		# Update the score hash
		foreach my $id (@$starts)
		{
			@ids_without = grep { $_ != $id } @$starts;
			$$score_hash{$id}{$_}{starts_close} = 1 for @ids_without;
		}

		my $ends = &ends_in_window($dbh, $chr, $left, $right);
		foreach my $id (@$ends)
		{
			@ids_without = grep { $_ != $id } @$ends;
			$$score_hash{$id}{$_}{ends_close} = 1 for @ids_without;
		}

		my $summits = &summits_in_window($dbh, $chr, $left, $right);
		foreach my $id (@$ends)
		{
			@ids_without = grep { $_ != $id } @$summits;
			$$score_hash{$id}{$_}{summits_close} = 1 for @ids_without;
		}
		
		$left += $increment;
		$right += $increment;
	}
}

sub starts_in_window {
	my ($dbh, $chr, $left, $right) = @_;
	my $query = "SELECT id FROM peaks WHERE chromosome=? AND start BETWEEN ? AND ?";
	my $sth = $dbh->prepare($query);
	$sth->execute($chr, $left, $right);
	
	my @starts;
	while( my $row = $sth->fetchrow_hashref())
	{
		push(@starts, $$row{id});
	}
	return \@starts;
}

sub ends_in_window {
	my ($dbh, $chr, $left, $right) = @_;
	my $query = "SELECT id FROM peaks WHERE chromosome=? AND stop BETWEEN ? AND ?";
	my $sth = $dbh->prepare($query);
	$sth->execute($chr, $left, $right);
	
	my @ends;
	while( my $row = $sth->fetchrow_hashref())
	{
		push(@ends, $$row{id});
	}
	return \@ends;
}

sub summits_in_window {
	my ($dbh, $chr, $left, $right) = @_;
	my $query = "SELECT id FROM peaks WHERE chromosome=? AND summit BETWEEN ? AND ?";
	my $sth = $dbh->prepare($query);
	$sth->execute($chr, $left, $right);
	
	my @summits;
	while( my $row = $sth->fetchrow_hashref())
	{
		push(@summits, $$row{id});
	}
	return \@summits;
}