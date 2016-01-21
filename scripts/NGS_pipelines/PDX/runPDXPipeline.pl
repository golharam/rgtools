#!/bin/env perl
use strict;

my $VERSION = "0.1";
my %SAMPLES;
my $dryRun = 0;

sub readSamples {
        my ($sampleSheet) = @_;
        print STDERR "Reading $sampleSheet...";

        open(IN, "<$sampleSheet") || die "Unable to open $sampleSheet\n";
        while (<IN>) {
                next if ($_ =~ m/^#/);
                chomp;
                my @fields = split(/\t/, $_);
                # Make sure file is properly formatted
                if (scalar(@fields) < 3) {
                        print STDERR "\nExpect 3 columns, but found ", scalar(@fields), "\n";
                        exit(-1);
                }
                # Make sure sample isn't already defined
                if (defined($SAMPLES{$fields[0]})) {
                        print STDERR "\n$fields[0] already defined.\n";
                        exit(-1);
                }
                $SAMPLES{$fields[0]}{'fq1'} = $fields[1];
                $SAMPLES{$fields[0]}{'fq2'} = $fields[2];
        }
        close(IN);
        print STDERR "\nFound ", scalar(keys %SAMPLES), " samples.\n";
}

sub runSamples {
        for my $sampleName (sort keys %SAMPLES) {
                next if (-d $sampleName);
                print STDERR "Submitting $sampleName...\n";

                my ($fq1, $fq2) = ($SAMPLES{$sampleName}{'fq1'}, $SAMPLES{$sampleName}{'fq2'});
                my $qsubCommand = "qsub -N $sampleName -v SAMPLE=$sampleName,FASTQ1=$fq1,FASTQ2=$fq2 pdx_pipeline.sh";
                if ($dryRun == 1) {
                        print "$qsubCommand\n";
                } else {
                        $_ = `$qsubCommand`;
                        $_ =~ m/Your job (\d+)/;
                        if (!defined($1)) {
                                print STDERR "\nUnable to determine job ID: $_";
                                exit(-1);
                        }
                        $SAMPLES{$sampleName}{'job'} = $1;
                        print STDERR "$1\n";
                }
        }
}

readSamples($ARGV[0]);
runSamples();
