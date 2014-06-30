### Try to generate a MOSAICS fit
### will vary the bgEst method and stuff if it fails the first time
sub try_fit
{
	my $mos = shift;
	eval { $mos->fit(); };
	if($@)
	{
		my $return = &vary_bgEst($mos);
		if($return) { return $return; }

		say "Could not generate the fit using bgEst values\nAttempting truncProb value changing";
		my @truncProb = (0.999, 0.995, 0.99, 0.85);
		for my $prob (@truncProb)
		{
			eval { $mos->fit({'truncProb' => $prob }); };
			if($@)
			{
				$mos->fit_name('');
			}
		}
	}
	if($mos->fit_name ne '') {
		return 1;
	} else {
		return -1;
	}
}

### Helper method for try_fit
sub vary_bgEst
{
	my $mos = shift;
	my @bgEst = qw|matchLow rMOM|;
	for my $opt (@bgEst)
	{
		eval { $mos->fit({'bgEst' => $opt}); };
		if($@)
		{
			say $@;
		} else {
			return $fit_name;
		}
	}
	return -1;
}