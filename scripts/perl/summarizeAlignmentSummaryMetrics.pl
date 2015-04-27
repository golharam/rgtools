#!/usr/bin/env perl
use strict;

my %SAMPLES;

# 1.  Read list of samples
open(IN, "<samples.txt") || die "Unable to open samples.txt\n";
while (<IN>) {
	chomp;
	my ($sample) = split(/\t/);

	getSampleAlignmentMetrics($sample);
}
close(IN);
printSummary();
exit(0);

sub getSampleAlignmentMetrics {
	my ($sample) = @_;

	open(ALN, "<analysis/qaqc/$sample.alignment_summary_metrics.txt") || die "Unable to open analysis/qaqc/$sample.alignment_summary_metrics.txt\n";
	my @metrics = <ALN>;
	close(ALN);
	chomp @metrics;

	my @fields = split(/\t/, $metrics[6]);
	my @values = split(/\t/, $metrics[12]);

	for (my $i = 0; $i < (scalar(@fields)); $i++) {
		$SAMPLES{$sample}{$fields[$i]} = $values[$i];
	}

}

sub printSummary {
	my @fields;

	# Print the header
	for my $sample (sort keys %SAMPLES) {
		print "\t$sample";
		@fields = sort keys $SAMPLES{$sample};
	}
	print "\n";

	for my $field (@fields) {
		print "$field";
		for my $sample (sort keys %SAMPLES) {
			print "\t".$SAMPLES{$sample}{$field};
		}
		print "\n";
	}
	print "\n";
}
