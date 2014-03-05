use Sequence;

my $seq = Sequence->new(
	length => 5,
	seq => "ATCGA",
);

print $seq->length;
my $ref = $seq->seq_array;
foreach my $nt (@$ref)
{
	#print "$nt\n";
}

#print $seq->base_at(1);
print $seq->range(1,3);