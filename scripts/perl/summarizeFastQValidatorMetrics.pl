#!/usr/bin/env perl
use strict;

my %SAMPLES;

# 1.  Read list of samples
open(IN, "<$ARGV[0]") || die "Unable to open $ARGV[0]\n";
while (<IN>) {
	chomp;
	my ($sample) = split(/\t/);

	getSampleFastqValidatorMetrics($sample);
}
close(IN);
printSummary();
exit(0);

sub getSampleFastqValidatorMetrics {
	my ($sample) = @_;

	open(F, "<analysis/qaqc/$sample.fq1.fastqValidator.txt") || die "Unable to open analysis/qaqc/$sample.fq1.fastqValidator.txt\n";
	my @metrics = <F>;
	close(F);
	chomp @metrics;

	$metrics[0] =~ m/with \d+ lines containing (\d+) sequences/;
	$SAMPLES{$sample}{'fq1_sequences'} = $1;

	$metrics[2] =~ m/Returning: \d+ : (\w+)/;
	$SAMPLES{$sample}{'fq1_status'} = $1;

	if (-e "analysis/qaqc/$sample.fq2.fastqValidator.txt") {
	        open(F, "<analysis/qaqc/$sample.fq2.fastqValidator.txt") || die "Unable to open analysis/qaqc/$sample.fq2.fastqValidator.txt\n";
	        my @metrics = <F>;
        	close(F);
        	chomp @metrics;

        	$metrics[0] =~ m/with \d+ lines containing (\d+) sequences/;
        	$SAMPLES{$sample}{'fq2_sequences'} = $1;

        	$metrics[2] =~ m/Returning: \d+ : (\w+)/;
        	$SAMPLES{$sample}{'fq2_status'} = $1;
	}
}

sub printSummary {
	my @fields;

	# Get the first sample to get a list of fields
	for my $sample (sort keys %SAMPLES) {
		@fields = sort keys $SAMPLES{$sample};
	}

	# Print the header
	for my $field (@fields) {
		print "\t$field";
	}
	print "\n";

	for my $sample (sort keys %SAMPLES) {
		print "$sample";
		for my $field (@fields) {
			print "\t".$SAMPLES{$sample}{$field};
		}
		print "\n";
	}
	print "\n";
}
