#!/usr/bin/env perl
use strict;

# This script creates the sample annotation file necessary for Xpress.  It uses
# the output of various other tools and collects them into one matrix.

my %SAMPLES;
my %ANNOTATION;

# 1.  Read list of samples
open(IN, "<samples.txt") || die "Unable to open samples.txt\n";
while (<IN>) {
	chomp;
	my ($sample) = split(/\t/);
	#$SAMPLES{$sample} = 1; # placeholder to record sample

	# todo: Raw read count
	
	# Samtools flagstat metrics
	$ANNOTATION{'mappedReads'} = "Covariate";
	collectSamtoolsFlagstatMetrics($sample);
}
close(IN);

# Output sample annotation
open(XPRESS, ">xpress/sample.annotation") || die "Unable to open xpress/sample.annotation\n";

my @factors = sort keys %ANNOTATION;

# Header Line 1
print XPRESS "X-Object File v";
for my $factor (@factors) {
	print XPRESS "\t$factor";
}
print XPRESS "\t<- Factor Name\n";

# Header Line 2
print XPRESS "Factor Type ->";
for my $factor (@factors) {
	print XPRESS "\t".$ANNOTATION{$factor};
}
print XPRESS "\tX-Object Name v\n";

for my $sample (sort keys %SAMPLES) {
	print XPRESS $sample;
	for my $factor (@factors) {
		print XPRESS "\t".$SAMPLES{$sample}{$factor};
	}
	print XPRESS "\t$sample\n";
}
close(XPRESS);

exit 0;

####
sub collectSamtoolsFlagstatMetrics {
	my ($sample) = @_;

	open(FLAGSTAT, "<$sample/$sample.genome.aln.bam.flagstat") ||
		print STDERR "Unable to open $sample/$sample.genome.aln.bam.flagstat\n";
	my @lines = <FLAGSTAT>;
	chomp @lines;
	close(FLAGSTAT);

	if ($lines[2] =~ /^(\d+)/) {
		$SAMPLES{$sample}{'mappedReads'} = $1;
	} else {
		$SAMPLES{$sample}{'mappedReads'} = -1;
		print STDERR "Unable to determine mapped read count for $sample\n";
	}
}
