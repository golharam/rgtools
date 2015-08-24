#!/usr/bin/env perl
use strict;

my %SAMPLES;

# 1.  Read list of samples
open(IN, "<$ARGV[0]") || die "Unable to open $ARGV[0]\n";
while (<IN>) {
	chomp;
	my ($sample) = split(/\t/);

	getSampleDedupMetrics($sample);
}
close(IN);
printSummary();
exit(0);

sub getSampleDedupMetrics {
	my ($sample) = @_;

	open(DEDUP, "<analysis/qaqc/$sample.dedup.metrics") || die "Unable to open analysis/qaqc/$sample.dedup.metrics\n";
	my @metrics = <DEDUP>;
	close(DEDUP);
	chomp @metrics;

	my @fields = split(/\t/, $metrics[6]);
	my @values = split(/\t/, $metrics[7]);

	for (my $i = 0; $i < (scalar(@fields)); $i++) {
		$SAMPLES{$sample}{$fields[$i]} = $values[$i];
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
