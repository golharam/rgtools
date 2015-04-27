#!/bin/env perl
use strict;
use warnings;

# Use this script to convert Illumina FASTQ headers from new style to old.
# Ex. gunzip -dc reads.fq.gz | fillIlluminaHeaders.pl > reads.fixed.fq

my $line = 0;
while (<STDIN>) {
	if ($line % 4 == 0) {
		chomp;
		my @fields = split(/\s/);
		#my $pair = ($fields[1] =~ m/^1:/) ? '1' : '2'
		my $pair;
		if ($fields[1] =~ m/^1:/) {
			$pair = '1';
		} elsif ($fields[1] =~ m/^[34]:/) { # not sure why this would be 3 or 4 instead of 2, but it is
			$pair = '2';
		} else {
			print STDERR "Unable to parse read name: $_\n";
			exit(-1);
		}
		print "$fields[0]".'#0/'."$pair\n";
	} else {
		print $_;
	}
	$line++;
}
