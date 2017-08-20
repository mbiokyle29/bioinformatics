#!/usr/bin/env perl
#
# Kyle McChesney
# filedo.t
# 
#
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use feature 'say';
use FileDo;
use Test::Simple tests => 4;

ok(my $test = FileDo->new(basedir => '/home/kyle/lab/FileDo/test/'), "creationg ok");
ok(my $test_reg = FileDo->new(basedir => '/home/kyle/lab/FileDo/test/', filetype => 'fasta'), "creating regex file okay");
ok($test_reg->print() == 0, "FileDo printing all the files");
ok($test_reg->do("cat (?)") == 0, "Manually cating the fasta file with the do method");
ok($test->do_all("cat (?)") == 0, "Do al");
