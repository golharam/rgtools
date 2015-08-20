#!/usr/bin/env perl
use strict;
use File::Basename;

# This script creates a smple sheet of (sample, fq1, fq2)
# Input is a text file of fq.gz.  FQ files must be named <sample>_[12].fq.gz

my %SAMPLES;

my $cwd = `pwd`;
chomp $cwd;

# allsamples.txt should be a text file of all the fastq.gz files
# Typically made with: find <path> -name '*.fastq.gz
for my $gzfile (`cat allsamples.txt`) {
	chomp $gzfile;

	# Breakup $gzfile into its parts
	my @suffixes = (".clipped.fastq.gz", ".fq.gz", ".fastq.gz");
	my ($filename, $dir, $suffix) = fileparse($gzfile, @suffixes);

	# Get the sample name from the filename
	if ($filename =~ /^(.*)_1$/) {
		my $sample = $1;
		if (defined($SAMPLES{$sample}{'fq1'})) {
			print STDERR "$sample fq1 already exists\n";
		} else {
			$SAMPLES{$sample}{'fq1'} = "$gzfile";
			$SAMPLES{$sample}{'dir'} = $dir;
		}
	} elsif ($filename =~ /^(.*)_2$/) {
		my $sample = $1;
		if (defined($SAMPLES{$sample}{'fq2'})) {
			print STDERR "$sample fq2 already exists\n";
		} else {
			$SAMPLES{$sample}{'fq2'} = "$gzfile";
		}
	}
}

for my $sample (keys %SAMPLES) {
	print join("\t", $sample, $SAMPLES{$sample}{'fq1'}, $SAMPLES{$sample}{'fq2'}), "\n";
}
