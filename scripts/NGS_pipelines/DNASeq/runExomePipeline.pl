#!/bin/env perl
use strict;
use warnings;

use File::Basename;

my $VERSION = "0.1";

my $dryRun = 0;
my %SAMPLES;

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

sub runStage1 {
	my ($project, $targetRegions) = @_;

	# Get the directory of where this script is.  
	my $dirname = dirname(__FILE__);

	for my $sampleName (sort keys %SAMPLES) {
		next if (-d $sampleName);
		print STDERR "Submitting $sampleName...\n";
		
		my ($fq1, $fq2) = ($SAMPLES{$sampleName}{'fq1'}, $SAMPLES{$sampleName}{'fq2'});
		my $qsubCommand = "qsub -N $sampleName -v SAMPLE=$sampleName,FASTQ1=$fq1,FASTQ2=$fq2,PROJECT=$project,AWS=1,TARGETREGIONS=$targetRegions $dirname/ExomePipeline.sh";
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

sub waitForSamples {
	for my $sampleName (sort keys %SAMPLES) {
		next if (!defined($SAMPLES{$sampleName}{'job'}));
		my $job = $SAMPLES{$sampleName}{'job'};
		print STDERR "Waiting for $sampleName ($job)...";

		my $state = getJobStatus($job);
		while ($state !~ m/c/) {
			sleep(60);
			$state = getJobStatus($job);
		}
		print STDERR "Done.\n";
	}
}

sub main {
	my $project = '';
	my $targetRegions = '';
	my $samplesFile = '';

	if (scalar(@ARGV) == 3) {
		$project = $ARGV[0];
		$targetRegions = $ARGV[1];
		$samplesFile = $ARGV[2];
	} else {
		print STDERR "Usage: $0 <project> <targetregions> <samples.txt>\n";
		exit(-1);
	}

	readSamples($samplesFile);
	runStage1($project, $targetRegions);
	waitForSamples();
}
main();

