package Query;
use Mouse;
use namespace::autoclean;
# A worth while abstraction

sub find_closest_peak_below {
	my ($self, $point, $canidates, $param) = @_;
	my $closest_peak = shift(@$canidates);
	my $closest_dist = $point - $closest_peak->$param;
	foreach my $cand_peak (@$canidates) {
		my $dist = $point - $cand_peak->$param;
		if($dist < $closest_dist && $dist >= 0) {

		}
	}
}




__PACKAGE__->meta->make_immutable;
1;