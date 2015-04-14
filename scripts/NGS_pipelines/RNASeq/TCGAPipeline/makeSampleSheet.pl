#!/usr/bin/env perl
use strict;

my %SAMPLES;

my $cwd = `pwd`;
chomp $cwd;

print "# Created by Ryan Golhar\n";
print "# file is used by my scripts to keep track of samples and file locations for pipeline\n";
print "# Basedir: $cwd";
print "#Sample\tFastQ1\tFastQ2\n";

for my $gzfile (`find /net/kraken/ng24/CA209-038_Biopsy_Pre_And_Post_Treatment_Plate_1_And_Plate_2_RNAseq -name '*.clipped.fastq.gz'`) {
	chomp $gzfile;
	my $baseName = `basename $gzfile`;
	$baseName =~ m/\w+-Plate-\d-(.+).clipped.fastq.gz/;
	my $sample = $1;
	if ($sample =~ /^(.*)_1$/) {
		my $sample_name = $1;
		if (defined($SAMPLES{$sample_name}{'fq1'})) {
			print STDERR "$sample_name fq1 already exists\n";
		} else {
			$SAMPLES{$sample_name}{'fq1'} = "$gzfile";
		}
	} elsif ($sample =~ /^(.*)_2$/) {
                my $sample_name = $1;
		if (defined($SAMPLES{$sample_name}{'fq2'})) {
			print STDERR "$sample_name fq2 already exists\n";
		} else {
			$SAMPLES{$sample_name}{'fq2'} = "$gzfile";
		}
	}
}

for my $sample (keys %SAMPLES) {
	print join("\t", $sample, $SAMPLES{$sample}{'fq1'}, $SAMPLES{$sample}{'fq2'}), "\n";
}
