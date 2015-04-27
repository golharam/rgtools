#!/bin/env perl
use strict;

open(BED, "<$ARGV[0]");

my $targets = 0;
my $bp = 0;
while (<BED>) {
	$targets++;
	chomp;
	my ($chr, $start, $end) = split(/\t/);
	#print "$chr\t$start\t$end\t".($end - $start)."\n";
	$bp += ($end - $start);
}

print "$targets targets spanning $bp bps\n";

close (BED);
