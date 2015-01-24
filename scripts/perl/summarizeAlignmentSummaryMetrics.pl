#!/usr/bin/env perl
use strict;

open(FILE, "<$ARGV[0]") || die "Unable to open $ARGV[0]";

my $printHeader = 0;
my @colNames;
my $sampleColIndex;
while (<FILE>) {
	next if length($_) == 0;
	chomp;

	# Collect the alignment metrics for all the samples
	if (m/^## METRICS CLASS/) {
		$_ = <FILE>;
		if ($printHeader == 0) {
			print $_;
			$printHeader = 1;
		}
		chomp;
		@colNames = split(/\t/);
		$sampleColIndex = 0;
		while ($colNames[$sampleColIndex] !~ m/SAMPLE/) { $sampleColIndex++ }

		my $line = <FILE>;		
		while (length($line) != 0) {
			chomp $line;
			my @fields = split(/\t/, $line);
			if (length($fields[$sampleColIndex]) > 0) {
				my $sample = $fields[$sampleColIndex];
				print "$line\n";		
			}
			$line = <FILE>;
		}
	}
}
close(FILE);
