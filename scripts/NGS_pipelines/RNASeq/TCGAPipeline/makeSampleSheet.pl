#!/usr/bin/env perl
use strict;

my %SAMPLES;

my $cwd = `pwd`;
chomp $cwd;

print STDERR "Using Source Directory: $ARGV[0]\n";
print STDERR "Using extension: $ARGV[1]\n";
my $srcdir = $ARGV[0];
my $ext = $ARGV[1];

print "# Created by Ryan Golhar\n";
print "# file is used by my scripts to keep track of samples and file locations for pipeline\n";
print "# Basedir: $ARGV[0]\n";
print "#Sample\tFastQ1\tFastQ2\n";

for my $file (`find $srcdir`) { 
	chomp $file;
	next if (index($file, $ext) == -1);

	my $baseName = `basename $file $ext`;
	#$baseName =~ m/\w+-Plate-\d-(.+).clipped.fastq.gz/;
	my $sample = $1;
	if ($sample =~ /^(.*)_1$/) {
		my $sample_name = $1;
		if (defined($SAMPLES{$sample_name}{'fq1'})) {
			print STDERR "$sample_name fq1 already exists\n";
		} else {
			$SAMPLES{$sample_name}{'fq1'} = "$file";
		}
	} elsif ($sample =~ /^(.*)_2$/) {
                my $sample_name = $1;
		if (defined($SAMPLES{$sample_name}{'fq2'})) {
			print STDERR "$sample_name fq2 already exists\n";
		} else {
			$SAMPLES{$sample_name}{'fq2'} = "$file";
		}
	}
}

for my $sample (keys %SAMPLES) {
	print join("\t", $sample, $SAMPLES{$sample}{'fq1'}, $SAMPLES{$sample}{'fq2'}), "\n";
}
