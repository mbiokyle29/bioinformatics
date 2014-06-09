#!/usr/bin/perl
use warnings;
use strict;
use Inline C => <<'END_C';

void greet() {
	printf("Hello World\n");
}
END_C
&greet;