#!/usr/bin/env perl
use strict;

my %SAMPLES;
my @FIELDS;

# 1.  Read list of samples
open(IN, "<$ARGV[0]") || die "Unable to open $ARGV[0]\n";
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

	my $i = 0;
	while ($i < scalar(@metrics)) {
		if ($metrics[$i] =~ m/^CATEGORY/) {
			@FIELDS = split(/\t/, $metrics[$i]);
		}
		$i++;
	}

	my @values;
	$i = scalar(@metrics) - 1;
	while ($i >= 0) {
		if (length($metrics[$i]) > 0) {
			@values = split(/\t/, $metrics[$i]);
			$i = -1;
		}
		$i--;
	}

	for (my $i = 0; $i < (scalar(@FIELDS)); $i++) {
		$SAMPLES{$sample}{$FIELDS[$i]} = $values[$i];
	}

}

sub printSummary {

	# Print the header
        for my $field (@FIELDS) {
                print "\t$field";
        }
        print "\n";

        for my $sample (sort keys %SAMPLES) {
                print "$sample";
                for my $field (@FIELDS) {
                        print "\t".$SAMPLES{$sample}{$field};
                }
                print "\n";
        }
        print "\n";
}
