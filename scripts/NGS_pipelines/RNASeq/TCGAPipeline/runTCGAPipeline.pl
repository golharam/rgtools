#!/bin/env perl
use strict;
use warnings;

use File::Basename;

my %SAMPLES;

sub readSamples {
	my ($sampleSheet) = @_;
	print STDERR "Reading $sampleSheet...";
	
	open(IN, "<$sampleSheet") || die "Unable to open $sampleSheet\n";
	while (<IN>) {
		chomp;
		my @fields = split(/\t/, $_);
		# Make sure file is properly formatted
		if (scalar(@fields) != 3) {
			print STDERR "\nExpect 3 columns, but found ", scalar(@fields), "\n";
			exit(-1);
		}
		my ($sampleName, $fq1, $fq2) = @fields;
		# Make sure sample isn't already defined
		if (defined($SAMPLES{$sampleName})) {
			print STDERR "\n$sampleName already defined.\n";
			exit(-1);
		}
		$SAMPLES{$sampleName}{'fq1'} = $fq1;
		$SAMPLES{$sampleName}{'fq2'} = $fq2;		
	}
	close(IN);
	print STDERR "\nFound ", scalar(keys %SAMPLES), " samples.\n";
}

sub runSamples {
	# Get the directory of where this script is.  
	my $dirname = dirname(__FILE__);

	for my $sampleName (sort keys %SAMPLES) {
		print STDERR "Submitting $sampleName...";
		
		my ($fq1, $fq2) = ($SAMPLES{$sampleName}{'fq1'}, $SAMPLES{$sampleName}{'fq2'});
		$_ = `qsub -N $sampleName -v SAMPLE=$sampleName,FASTQ1=$fq1,FASTQ2=$fq2 $dirname/unc_rnaseqV2_pipeline.v2.sh`;
		$_ =~ m/Your job (\d+)/;
		if (!defined($1)) {
			print STDERR "\nUnable to determine job ID: $_";
			exit(-1);
		}
		$SAMPLES{$sampleName}{'job'} = $1;
		print STDERR "$1\n";
	}
}

sub getJobStatus {
	my ($job) = @_;
	$_ = `qstat | grep $job`;
	chomp;
	$_ =~ s/^\s+|\s+$//g;	# remove leading and trailing spaces
	my @fields = split(/\s+/, $_);
	my $state = $fields[4];

	return 'c' if (length($_) == 0);
	return $state;
}

sub waitForSamples() {
	for my $sampleName (sort keys %SAMPLES) {
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

sub makeExpressionMatrices {
	# Get the directory of where this script is.  
	my $dirname = dirname(__FILE__);

	# Gene-level expression
	my $abundanceFiles = "";
	for my $sampleName (sort keys %SAMPLES) {
		my $abundanceFile = "$sampleName/$sampleName.rsem.genes.results";
		if (! -e $abundanceFile) {
			die "Unable to locate $abundanceFile\n";
		}
		$abundanceFiles .= "$abundanceFile ";
	}
	`$dirname/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix genes $abundanceFiles`;

	# Isoform-level expression
	$abundanceFiles = "";
	for my $sampleName (sort keys %SAMPLES) {
		my $abundanceFile = "$sampleName/$sampleName.rsem.isoforms.results";
		if (! -e $abundanceFile) {
			die "Unable to locate $abundanceFile\n";
		}
		$abundanceFiles .= "$abundanceFile ";
	}
	`$dirname/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix isoforms $abundanceFiles`;

	`$dirname/TCGAQuartileNormalizationUNC.pl -s 1 -c -1 -o genes.UQ.matrix genes.counts.matrix`;
	`$dirname/TCGAQuartileNormalizationUNC.pl -s 1 -c -1 -o isoforms.UQ.matrix isoforms.counts.matrix`;
}

sub main {
	if (scalar(@ARGV) != 1) {
		print STDERR "Usage: $0 <samples.txt>\n";
		exit(-1);
	}

	readSamples($ARGV[0]);
	runSamples();
	waitForSamples();
	makeExpressionMatrices();
}

main();
