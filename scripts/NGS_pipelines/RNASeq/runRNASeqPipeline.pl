#!/bin/env perl
use strict;
use warnings;

=head1 SYNOPSIS

runRNASeqQCPipeline.pl
	[--dryRun                (default: False] 
        [--outdir <output dir>   (default: analysis/{samples})
	[--refmodel <refmodel>   (default: hg19ERCC.refSeq)]
	[--tmpdir <tmpdir>       (default: /scratch)] 
	--samples <samples.txt>

=cut

use File::Basename;
use Getopt::Long;
use Pod::Usage 'pod2usage';

my $VERSION = "0.1";

my $dryRun = 0;
my %SAMPLES;

sub main {
	my $outdir = `pwd`; chomp $outdir;
	my $refmodel = 'hg19ERCC.refSeq';
	my $tmpdir = '/scratch';
	my $sampleSheet = '';

	# TBD: Something is wrong with 
	# 'tmpdir:s' => $tmpdir
	GetOptions('dryRun' => \$dryRun,
		   'outdir=s' => \$outdir,
		   'refmodel=s' => \$refmodel,
		   'tmpdir=s' => \$tmpdir,
		   'samples=s' => \$sampleSheet);
	pod2usage unless $sampleSheet;

	print STDERR "dryRun: $dryRun\n";
	print STDERR "outdir: $outdir\n";
	print STDERR "refmodel: $refmodel\n";
	print STDERR "tpmdir: $tmpdir\n";
	print STDERR "samples: $sampleSheet\n";
	print STDERR "\n";

	readSamples($sampleSheet);
	runSamples($refmodel, $tmpdir, $outdir);
	if ($dryRun != 1) {
		waitForSamples();
		makeExpressionMatrices($outdir);
	}
}

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
	my ($refmodel, $tmpdir, $outdir) = @_;

	for my $sampleName (sort keys %SAMPLES) {
		next if (-d "$outdir/$sampleName");
		
		my ($fq1, $fq2) = ($SAMPLES{$sampleName}{'fq1'}, $SAMPLES{$sampleName}{'fq2'});
		my $qsubCommand = "qsub -N $sampleName -pe orte 8 -v SAMPLE=$sampleName,FASTQ1=$fq1,FASTQ2=$fq2,TMP_DIR=$tmpdir,REFMODEL=$refmodel,OUTDIR=$outdir $dirname/rnaseq_pipeline.sh";
		if ($dryRun == 1) {
			print "$qsubCommand\n";
		} else {
			print STDERR "Submitting $sampleName...";
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
	my ($outdir) = @_;

	# Gene-level expression
	my $abundanceFiles = "";
	for my $sampleName (sort keys %SAMPLES) {
		my $abundanceFile = "$outdir/$sampleName/$sampleName.genes.results";
		if (! -e $abundanceFile) {
			die "Unable to locate $abundanceFile\n";
		}
		$abundanceFiles .= "$abundanceFile ";
	}
	`$dirname/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix $outdir/genes $abundanceFiles`;

	# Isoform-level expression
	$abundanceFiles = "";
	for my $sampleName (sort keys %SAMPLES) {
		my $abundanceFile = "$outdir/$sampleName/$sampleName.isoforms.results";
		if (! -e $abundanceFile) {
			die "Unable to locate $abundanceFile\n";
		}
		$abundanceFiles .= "$abundanceFile ";
	}
	`$dirname/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix $outdir/isoforms $abundanceFiles`;

	`$dirname/TCGAQuartileNormalizationUNC.pl -s 1 -c -1 -o $outdir/genes.UQ.matrix $outdir/genes.counts.matrix`;
	`$dirname/TCGAQuartileNormalizationUNC.pl -s 1 -c -1 -o $outdir/isoforms.UQ.matrix $outdir/isoforms.counts.matrix`;

	# TBD: Make samples manifest for Xpress
}

main();