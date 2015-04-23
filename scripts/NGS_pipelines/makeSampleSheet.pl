#!/usr/bin/env perl
use strict;
use File::Basename;

# This script creates a smple sheet of (sample, fq1, fq2)
# Input is a text file of fq.gz.  FQ files must be named <sample>_[12].fq.gz

my %SAMPLES;

my $cwd = `pwd`;
chomp $cwd;

for my $gzfile (`find F14FTSUSAT1059/ -name '*.fq.gz'`) {
	chomp $gzfile;

	# Breakup $gzfile into its parts
	my @suffixes = (".fq.gz", ".fastq.gz");
	my ($filename, $dir, $suffix) = fileparse($gzfile, @suffixes);

	# Get the sample name from the filename
	if ($filename =~ /^(\w+)_pe_1$/) {
		my $sample = $1;
		if (defined($SAMPLES{$sample}{'fq1'})) {
			print STDERR "$sample fq1 already exists\n";
		} else {
			$SAMPLES{$sample}{'fq1'} = "$cwd/$gzfile";
		}
	} elsif ($filename =~ /^(.*)_pe_2$/) {
                my $sample = $1;
		if (defined($SAMPLES{$sample}{'fq2'})) {
			print STDERR "$sample fq2 already exists\n";
		} else {
			$SAMPLES{$sample}{'fq2'} = "$cwd/$gzfile";
		}
	}
}

for my $sample (keys %SAMPLES) {
	print join("\t", $sample, $SAMPLES{$sample}{'fq1'}, $SAMPLES{$sample}{'fq2'}), "\n";
}
