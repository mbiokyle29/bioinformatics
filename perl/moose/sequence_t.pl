use Sequence;

my $seq = Sequence->new(
	seq => "ATCGA"
);
print $seq->seq;
print $seq->length."\n";
my $ref = $seq->seq_arr;

foreach my $nt (@$ref)
{
	print "$nt";
}
print "\n";

print "Range 1 to 2: \n";
print $seq->range(1,2); 
print $seq->base_at($seq->length);
