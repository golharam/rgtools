#!/bin/env perl
use strict;
use warnings;

use File::Basename;

my %SAMPLES;
# Order of fields by line:
# FastQValidator
# FastQC
my @FASTQVALIDATOR_METRICS = ('fq1RawCount' , 'fq2RawCount');
my @FASTQC_METRICS = ('sequenceLength', 'pctGC');
my @CONTAMINATION_METRICS = ('contaminatedReadPairs', 'uncontaminatedReadPairs');
my @PICARD_ALNMETRICS = ('TOTAL_READS', 'PF_READS', 'PF_NOISE_READS', 'PF_READS_ALIGNED', 'PF_ALIGNED_BASES', 'PF_HQ_ALIGNED_READS', 'PF_HQ_ALIGNED_BASES', 'PF_HQ_ALIGNED_Q20_BASES', 'PF_HQ_MEDIAN_MISMATCHES', 'PF_MISMATCH_RATE', 'PF_HQ_ERROR_RATE', 'PF_INDEL_RATE', 'MEAN_READ_LENGTH', 'READS_ALIGNED_IN_PAIRS', 'PCT_READS_ALIGNED_IN_PAIRS', 'BAD_CYCLES', 'STRAND_BALANCE', 'PCT_CHIMERAS', 'PCT_ADAPTER');
my @PICARD_INSERTSIZEMETRICS = ('MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE', 'MAX_INSERT_SIZE', 'MEAN_INSERT_SIZE', 'STANDARD_DEVIATION', 'READ_PAIRS', 'PAIR_ORIENTATION', 'WIDTH_OF_10_PERCENT', 'WIDTH_OF_20_PERCENT', 'WIDTH_OF_30_PERCENT', 'WIDTH_OF_40_PERCENT', 'WIDTH_OF_50_PERCENT', 'WIDTH_OF_60_PERCENT', 'WIDTH_OF_70_PERCENT', 'WIDTH_OF_80_PERCENT', 'WIDTH_OF_90_PERCENT', 'WIDTH_OF_99_PERCENT');
my @PICARD_RNASEQMETRICS = ('PF_BASES', 'PF_ALIGNED_BASES', 'RIBOSOMAL_BASES', 'CODING_BASES', 'UTR_BASES', 'INTRONIC_BASES', 'INTERGENIC_BASES', 'IGNORED_READS', 'CORRECT_STRAND_READS', 'INCORRECT_STRAND_READS', 'PCT_RIBOSOMAL_BASES', 'PCT_CODING_BASES', 'PCT_UTR_BASES', 'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES', 'PCT_MRNA_BASES', 'PCT_USABLE_BASES', 'PCT_CORRECT_STRAND_READS', 'MEDIAN_CV_COVERAGE', 'MEDIAN_5PRIME_BIAS', 'MEDIAN_3PRIME_BIAS', 'MEDIAN_5PRIME_TO_3PRIME_BIAS');

sub readSamples {
	my ($sampleSheet) = @_;
	print STDERR "Reading $sampleSheet...";
	
	open(IN, "<$sampleSheet") || die "Unable to open $sampleSheet\n";
	while (<IN>) {
		next if ($_ =~ m/^#/);
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
		next if (-d $sampleName);
		print STDERR "Submitting $sampleName...";
		
		my ($fq1, $fq2) = ($SAMPLES{$sampleName}{'fq1'}, $SAMPLES{$sampleName}{'fq2'});
		$_ = `qsub -N $sampleName -v SAMPLE=$sampleName,FASTQ1=$fq1,FASTQ2=$fq2,SUBSAMPLE=500000 $dirname/rnaseqqc_pipeline.sh`;
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

sub collectMetrics {
	for my $sampleName (sort keys %SAMPLES) {
		my @data;

		# Collect FastQValidator Metrics
		open(F, "<$sampleName/$sampleName.fq1Validator.txt") || die "Unable to open $sampleName/$sampleName.fq1Validator.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		if ($data[0] =~ m/containing (\d+) sequences./) {
			$SAMPLES{$sampleName}{'fastQValidator'}{'fq1RawCount'} = $1;
		}
		if (-e "$sampleName/$sampleName.fq2Validator.txt") {
	                open(F, "<$sampleName/$sampleName.fq2Validator.txt") || die "Unable to open $sampleName/$sampleName.fq2Validator.txt";
	                @data = <F>;
        	        close(F);
                	chomp(@data);
	                if ($data[0] =~ m/containing (\d+) sequences./) {
        	                $SAMPLES{$sampleName}{'fastQValidator'}{'fq2RawCount'} = $1;
                	}
		}

		# Collect FastQC Metrics
		my $fastq_dir = $sampleName.'_1.fq_fastqc';
		if (-d "$sampleName/$fastq_dir") {
			open(F, "<$sampleName/$fastq_dir/fastqc_data.txt") || die "Unable to open $sampleName/$fastq_dir/fastqc_data.txt";
			@data = <F>;
			close(F);
			chomp(@data);

			if ($data[8] =~ m/Sequence length\t(\d+)/) {
				# TBD: This should come from picard::CollectAlignmentMetrics
				$SAMPLES{$sampleName}{'FastQC'}{'sequenceLength'} = $1;
			}
			if ($data[9] =~ m/\%GC\t(\d+)/) {
				$SAMPLES{$sampleName}{'FastQC'}{'pctGC'} = $1;
			}
			
			open(F, "<$sampleName/$fastq_dir/summary.txt") || die "Unable to open $sampleName/$fastq_dir/summary.txt";
			while (<F>) {
				chomp;
				my @data = split(/\t/);
				$SAMPLES{$sampleName}{'FastQC'}{$data[1]} = $data[0];
			}
			close(F);
		}

		# Collect Contamination Metrics
		open(F, "<$sampleName/contamination/$sampleName.contamination.log")
			|| die "Unable to open $sampleName/contamination/$sampleName.contamination.log\n";
		@data = <F>;
		close(F);
		chomp @data;
		my $totalReadPairs = $1 if ($data[0] =~ m/^(\d+)/);
		my $uncontaminatedReadPairs = $1 if ($data[2] =~ m/^\s+(\d+)/);
		my $contaminatedReadPairs = $totalReadPairs - $uncontaminatedReadPairs;
		$SAMPLES{$sampleName}{'contamination::contaminatedReadPairs'} = $contaminatedReadPairs;
		$SAMPLES{$sampleName}{'contamination::uncontaminatedReadPairs'} = $uncontaminatedReadPairs;

		# Collect Picard Alignment Metrics
		open(F, "<$sampleName/$sampleName.alnMetrics.txt") || die "Unable to open $sampleName/$sampleName.alnMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		my @fieldNames = split(/\t/, $data[6]);
		my @fieldValues = split(/\t/, $data[9]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'picard::AlignmentMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}

		# Collect Picard Insert Size Metrics
		open(F, "<$sampleName/$sampleName.insertSizeMetrics.txt") || die "Unable to open $sampleName/$sampleName.insertSizeMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		@fieldNames = split(/\t/, $data[6]);
		@fieldValues = split(/\t/, $data[7]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'picard::InsertSizeMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}

		# Collect Picard RNA-Seq Metrics
		open(F, "<$sampleName/$sampleName.rnaseqMetrics.txt") || die "Unable to open $sampleName/$sampleName.rnaseqMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		@fieldNames = split(/\t/, $data[6]);
		@fieldValues = split(/\t/, $data[7]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'picard::RNASeqMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}
	}
}

sub printMetrics {
	open(QC, "<qcMetrics.txt") || die "Unable to open qcMetrics.txt";
	
	for my $fields (@FASTQVALIDATOR_METRICS, @FASTQC_METRICS, @CONTAMINATION_METRICS, @PICARD_ALNMETRICS, @PICARD_INSERTSIZEMETRICS, @PICARD_RNASEQMETRICS) {
		print QC "\t$fields";
	}
	print QC "\n";

	for my $sampleName (sort keys %SAMPLES) {
		print QC "$sampleName";

		# Print out FastQValidator
		for my $fields (@FASTQVALIDATOR_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{'fastQValidator'}{$fields};
		}

		# Print out FastQC
		for my $fields (@FASTQC_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{'FastQC'}{$fields};
		}

		# Print out contamination metrics
		for my $fields (@CONTAMINATION_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{"contamination::$fields"};
		}
		
		# Print out Picard Alignment Metrics
		for my $fields (@PICARD_ALNMETRICS) {
			print QC "\t";
			if (defined($SAMPLES{$sampleName}{'picard::AlignmentMetrics'}{$fields})) {
				print QC $SAMPLES{$sampleName}{'picard::AlignmentMetrics'}{$fields};
			}
		}

		# Print out Picard Insert Size Metrics
		for my $fields (@PICARD_INSERTSIZEMETRICS) {
			print QC "\t";
			if (defined($SAMPLES{$sampleName}{'picard::InsertSizeMetrics'}{$fields})) {
				print QC $SAMPLES{$sampleName}{'picard::InsertSizeMetrics'}{$fields};
			}
		}

		# Print out Picard RNA-Seq Metrics
		for my $fields (@PICARD_RNASEQMETRICS) {
			print QC "\t";
			if (defined($SAMPLES{$sampleName}{'picard::RNASeqMetrics'}{$fields})) {
				print QC $SAMPLES{$sampleName}{'picard::RNASeqMetrics'}{$fields};
			}
		}

		print QC "\n";
	}
	print QC "\n";
	close(QC);
}

sub main {
	if (scalar(@ARGV) != 1) {
		print STDERR "Usage: $0 <samples.txt>\n";
		exit(-1);
	}

	readSamples($ARGV[0]);
	runSamples();
	waitForSamples();
	collectMetrics();
	printMetrics();
}

main();
