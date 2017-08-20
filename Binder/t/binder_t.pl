#!/usr/bin/perl
use warnings;
use strict;
use lib "/home/kyle/lab/Binder/lib/";
use Binder;
use Test::Simple tests => 10;
use Data::Printer;

# Basic one param
ok(my $binder = Binder->new(base_string => "This is a (?) String"), 'created object right');
ok($binder->base_string eq "This is a (?) String", 'Base string is right');
ok($binder->bind("TEST"), 'Binding works');
ok($binder->bound_string eq "This is a TEST String", "We bound!");
p($binder);

# Three params
ok($binder = Binder->new(base_string => "This is a (?) String with (?) (?)"), 'created second object right');
ok($binder->bind("TEST", "THREE", "PARAMS"), "Bound with three");
ok($binder->bound_string eq "This is a TEST String with THREE PARAMS", "Bound with three params");
p($binder);

# Three params + custom matcher
ok($binder = Binder->new(
	base_string => "This is a * String with * *",
	bind_symbol => "*"
), 'created second object right');
ok($binder->bind("TEST", "THREE", "PARAMS"), "Bound with three");
ok($binder->bound_string eq "This is a TEST String with THREE PARAMS", "Bound with three params and custom symbol");
p($binder);
