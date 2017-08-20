Binder
======

DBI like param binding on normal strings
````perl
use Binder;

# create a binder object, with a base string to be substitued
# the default binding symbol is (?)
my $binder = Binder->new(base_string => "This is a test string with a binder value: (?)");

# bound string
$binder->bind("NEW VALUE");
my $bound = $binder->bound_string();

# check slot number
if($binder->slot_number() == 1) {
  $binder->bind("another New value");
}

````
