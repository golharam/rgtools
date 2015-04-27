#!/bin/env perl
use strict;

my $targets = 0;
my $bp = 0;
while (<STDIN>) {
	$targets++;
	chomp;
	my ($chr, $start, $end) = split(/\t/);
	#print "$chr\t$start\t$end\t".($end - $start)."\n";
	$bp += ($end - $start);
}

print "$targets targets spanning $bp bps\n";
