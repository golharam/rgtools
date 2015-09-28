#!/bin/env perl
use strict;
use warnings;

use File::Basename;

my $VERSION = "0.5.3f";

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

sub runSamples {
	# Get the directory of where this script is.  
	my $dirname = dirname(__FILE__);
	my ($subsample, $species, $tmpdir) = @_;

	for my $sampleName (sort keys %SAMPLES) {
		next if (-d $sampleName);
		print STDERR "Submitting $sampleName...\n";
		
		my ($fq1, $fq2) = ($SAMPLES{$sampleName}{'fq1'}, $SAMPLES{$sampleName}{'fq2'});
		my $qsubCommand = "qsub -N $sampleName -v SAMPLE=$sampleName,FASTQ1=$fq1,FASTQ2=$fq2,USE_STAR=0,AWS=0,SUBSAMPLE=$subsample,TMP_DIR=$tmpdir,REFERENCE=$species $dirname/rnaseqqc_pipeline.sh";
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

sub collectAndPrintMetrics {
	my @FASTQVALIDATOR_METRICS = ('fq1RawCount' , 'fq2RawCount');
	my @FASTQC_METRICS = ('Basic Statistics', 'Per base sequence quality', 'Per tile sequence quality', 'Per sequence quality scores', 'Per base sequence content', 'Per sequence GC content', 'Per base N content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences', 'Adapter Content', 'Kmer Content', 'pctA', 'pctC', 'pctG', 'pctT', 'pctGC');
	my @CONTAMINATION_HEADER = ('ContaminatedReadPairs', '%Contamination');
	my @CONTAMINATION_METRICS = ('READS_ALIGNED_IN_PAIRS', '%Contamination');
	my @PICARD_ALNMETRICS = ('MEAN_READ_LENGTH', 'PF_NOISE_READS', 'PF_READS_ALIGNED', 'PF_ALIGNED_BASES', 'PF_HQ_ALIGNED_READS', 'PF_HQ_ALIGNED_BASES', 'PF_HQ_ALIGNED_Q20_BASES', 'PF_HQ_MEDIAN_MISMATCHES', 'PF_MISMATCH_RATE', 'PF_HQ_ERROR_RATE', 'PF_INDEL_RATE', 'READS_ALIGNED_IN_PAIRS', 'PCT_READS_ALIGNED_IN_PAIRS', 'BAD_CYCLES', 'STRAND_BALANCE', 'PCT_CHIMERAS', 'PCT_ADAPTER');
	my @PICARD_INSERTSIZEMETRICS = ('MEDIAN_INSERT_SIZE', 'MEDIAN_ABSOLUTE_DEVIATION', 'MIN_INSERT_SIZE', 'MAX_INSERT_SIZE', 'MEAN_INSERT_SIZE', 'STANDARD_DEVIATION', 'READ_PAIRS', 'PAIR_ORIENTATION');
	my @PICARD_RNASEQMETRICS = ('PF_BASES', 'PF_ALIGNED_BASES', 'RIBOSOMAL_BASES', 'CODING_BASES', 'UTR_BASES', 'INTRONIC_BASES', 'INTERGENIC_BASES', 'IGNORED_READS', 'CORRECT_STRAND_READS', 'INCORRECT_STRAND_READS', 'PCT_RIBOSOMAL_BASES', 'PCT_CODING_BASES', 'PCT_UTR_BASES', 'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES', 'PCT_MRNA_BASES', 'PCT_USABLE_BASES', 'PCT_CORRECT_STRAND_READS', 'MEDIAN_CV_COVERAGE', 'MEDIAN_5PRIME_BIAS', 'MEDIAN_3PRIME_BIAS', 'MEDIAN_5PRIME_TO_3PRIME_BIAS');
	my @HEMOGLOBIN_METRICS = ('HemogloginMappedReads');
	my @ERCC_METRICS = ('ERCC-00002','ERCC-00003','ERCC-00004','ERCC-00009','ERCC-00012','ERCC-00013','ERCC-00014','ERCC-00016','ERCC-00017','ERCC-00019','ERCC-00022','ERCC-00024','ERCC-00025','ERCC-00028','ERCC-00031',
			    'ERCC-00033','ERCC-00034','ERCC-00035','ERCC-00039','ERCC-00040','ERCC-00041','ERCC-00042','ERCC-00043','ERCC-00044','ERCC-00046','ERCC-00048','ERCC-00051','ERCC-00053','ERCC-00054','ERCC-00057',
			    'ERCC-00058','ERCC-00059','ERCC-00060','ERCC-00061','ERCC-00062','ERCC-00067','ERCC-00069','ERCC-00071','ERCC-00073','ERCC-00074','ERCC-00075','ERCC-00076','ERCC-00077','ERCC-00078','ERCC-00079',
			    'ERCC-00081','ERCC-00083','ERCC-00084','ERCC-00085','ERCC-00086','ERCC-00092','ERCC-00095','ERCC-00096','ERCC-00097','ERCC-00098','ERCC-00099','ERCC-00104','ERCC-00108','ERCC-00109','ERCC-00111',
			    'ERCC-00112','ERCC-00113','ERCC-00116','ERCC-00117','ERCC-00120','ERCC-00123','ERCC-00126','ERCC-00130','ERCC-00131','ERCC-00134','ERCC-00136','ERCC-00137','ERCC-00138','ERCC-00142','ERCC-00143',
			    'ERCC-00144','ERCC-00145','ERCC-00147','ERCC-00148','ERCC-00150','ERCC-00154','ERCC-00156','ERCC-00157','ERCC-00158','ERCC-00160','ERCC-00162','ERCC-00163','ERCC-00164','ERCC-00165','ERCC-00168',
			    'ERCC-00170', 'ERCC-00171');

	# Open output file
	open(QC, ">qcMetrics.txt") || die "Unable to open qcMetrics.txt";

	# Print header line
	for my $fields (@FASTQVALIDATOR_METRICS, @FASTQC_METRICS, @CONTAMINATION_HEADER, 
			@PICARD_ALNMETRICS, @PICARD_INSERTSIZEMETRICS, @PICARD_RNASEQMETRICS, 
			@HEMOGLOBIN_METRICS, @ERCC_METRICS) {
		print QC "\t$fields";
	}
	print QC "\n";

	for my $sampleName (sort keys %SAMPLES) {
		my @data;

		# Print Sample Name
		print QC "$sampleName";

		# Collect FastQValidator Metrics
		open(F, "<analysis/$sampleName/$sampleName.fq1Validator.txt") || die "Unable to open analysis/$sampleName/$sampleName.fq1Validator.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		if ($data[0] =~ m/containing (\d+) sequences./) {
			$SAMPLES{$sampleName}{'fastQValidator'}{'fq1RawCount'} = $1;
		}
		if (-e "analysis/$sampleName/$sampleName.fq2Validator.txt") {
	                open(F, "<analysis/$sampleName/$sampleName.fq2Validator.txt") || die "Unable to open analysis/$sampleName/$sampleName.fq2Validator.txt";
	                @data = <F>;
        	        close(F);
                	chomp(@data);
	                if ($data[0] =~ m/containing (\d+) sequences./) {
        	                $SAMPLES{$sampleName}{'fastQValidator'}{'fq2RawCount'} = $1;
                	}
		}
		# Print out FastQValidator
		for my $fields (@FASTQVALIDATOR_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{'fastQValidator'}{$fields};
		}

		# Collect FastQC Metrics
		my $fastq_dir = $sampleName.'_1.subsampled_fastqc';
		if (-d "analysis/$sampleName/$fastq_dir") {
			open(F, "<analysis/$sampleName/$fastq_dir/fastqc_data.txt") || die "Unable to open analysis/$sampleName/$fastq_dir/fastqc_data.txt";
			@data = <F>;
			close(F);
			chomp(@data);

			if ($data[8] =~ m/Sequence length\t(\d+)/) {
				# TBD: This should come from picard::CollectAlignmentMetrics
				$SAMPLES{$sampleName}{'FastQC'}{'sequenceLength'} = $1;
			}
			if ($data[9] =~ m/\%A\t(\d+)/) {
				$SAMPLES{$sampleName}{'FastQC'}{'pctA'} = $1;
			}
			if ($data[10] =~ m/\%C\t(\d+)/) {
				$SAMPLES{$sampleName}{'FastQC'}{'pctC'} = $1;
			}
			if ($data[11] =~ m/\%G\t(\d+)/) {
				$SAMPLES{$sampleName}{'FastQC'}{'pctG'} = $1;
			}
			if ($data[12] =~ m/\%T\t(\d+)/) {
				$SAMPLES{$sampleName}{'FastQC'}{'pctT'} = $1;
			}
			if ($data[13] =~ m/\%GC\t(\d+)/) {
				$SAMPLES{$sampleName}{'FastQC'}{'pctGC'} = $1;
			}
			
			open(F, "<analysis/$sampleName/$fastq_dir/summary.txt") || die "Unable to open analysis/$sampleName/$fastq_dir/summary.txt";
			while (<F>) {
				chomp;
				my @data = split(/\t/);
				$SAMPLES{$sampleName}{'FastQC'}{$data[1]} = $data[0];
			}
			close(F);
		}
		# Print out FastQC
		for my $fields (@FASTQC_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{'FastQC'}{$fields};
		}

		# Collect Contamination Metrics
		open(F, "<analysis/$sampleName/contamination/$sampleName.contaminated.alnMetrics.txt")
			|| die "Unable to open analysis/$sampleName/contamination/$sampleName.contaminated.alnMetrics.txt\n";
		@data = <F>;
		close(F);
		chomp @data;
		my @fieldNames = split(/\t/, $data[6]);
		my @fieldValues = split(/\t/, $data[7]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'contamination::AlignmentMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}
		$SAMPLES{$sampleName}{'contamination::AlignmentMetrics'}{'%Contamination'} = ($SAMPLES{$sampleName}{'contamination::AlignmentMetrics'}{'READS_ALIGNED_IN_PAIRS'} / $SAMPLES{$sampleName}{'fastQValidator'}{'fq1RawCount'})*100;
		# Print out contamination metrics
		for my $fields (@CONTAMINATION_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{"contamination::AlignmentMetrics"}{$fields};
		}

		# Collect Picard Alignment Metrics
		open(F, "<analysis/$sampleName/$sampleName.alnMetrics.txt") || die "Unable to open analysis/$sampleName/$sampleName.alnMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		@fieldNames = split(/\t/, $data[6]);
		@fieldValues = split(/\t/, $data[9]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'picard::AlignmentMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}
		# Print out Picard Alignment Metrics
		for my $fields (@PICARD_ALNMETRICS) {
			print QC "\t";
			if (defined($SAMPLES{$sampleName}{'picard::AlignmentMetrics'}{$fields})) {
				print QC $SAMPLES{$sampleName}{'picard::AlignmentMetrics'}{$fields};
			}
		}

		# Collect Picard Insert Size Metrics
		open(F, "<analysis/$sampleName/$sampleName.insertSizeMetrics.txt") || die "Unable to open analysis/$sampleName/$sampleName.insertSizeMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		@fieldNames = split(/\t/, $data[6]);
		@fieldValues = split(/\t/, $data[7]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'picard::InsertSizeMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}
		# Print out Picard Insert Size Metrics
		for my $fields (@PICARD_INSERTSIZEMETRICS) {
			print QC "\t";
			if (defined($SAMPLES{$sampleName}{'picard::InsertSizeMetrics'}{$fields})) {
				print QC $SAMPLES{$sampleName}{'picard::InsertSizeMetrics'}{$fields};
			}
		}

		# Collect Picard RNA-Seq Metrics
		open(F, "<analysis/$sampleName/$sampleName.rnaseqMetrics.txt") || die "Unable to open analysis/$sampleName/$sampleName.rnaseqMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		@fieldNames = split(/\t/, $data[6]);
		@fieldValues = split(/\t/, $data[7]);
		for (my $i = 0; $i < scalar(@fieldNames); $i++) {
			$SAMPLES{$sampleName}{'picard::RNASeqMetrics'}{$fieldNames[$i]} = $fieldValues[$i];
		}
		# Print out Picard RNA-Seq Metrics
		for my $fields (@PICARD_RNASEQMETRICS) {
			print QC "\t";
			if (defined($SAMPLES{$sampleName}{'picard::RNASeqMetrics'}{$fields})) {
				print QC $SAMPLES{$sampleName}{'picard::RNASeqMetrics'}{$fields};
			}
		}

		# Collect Hemoglobin Metrics
		open(F, "<analysis/$sampleName/$sampleName.hemoglobinMetrics.txt") || die "Unable to open analysis/$sampleName/$sampleName.hemoglobinMetrics.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		$data[4] =~ m/(\d+)$/;
		$SAMPLES{$sampleName}{'hemoglobin'}{'HemogloginMappedReads'} = $1;
		print QC "\t$1";
		
		# Collect ERCC Metrics
		# $SAMPLE.erccMetrics.txt has the format:
		# chr
		# start
		# end
		# The number of reads that overlapped an ERCC sequence.
		# The number of bases that had non-zero coverage.
		# The length of the ERCC entry.
		# The fraction of bases that had non-zero coverage.
		open(F, "<analysis/$sampleName/$sampleName.idxStats.txt") || die "Unable to open analysis/$sampleName/$sampleName.idxStats.txt";
		@data = <F>;
		close(F);
		chomp(@data);
		for (my $i = 0; $i < scalar(@data); $i++) {
			my @fields = split(/\t/, $data[$i]);	# chr, chr length, mapped read count, unmapped read cout
			if ($fields[0] =~ m/(ERCC-\d+)/) {
				$SAMPLES{$sampleName}{'erccReadCount'}{$1} = $fields[2];
			}
		}
		# Print out ERCC Metrics
		for my $fields (@ERCC_METRICS) {
			print QC "\t".$SAMPLES{$sampleName}{'erccReadCount'}{$fields};
		
		}

		print QC "\n";
	}
	close(QC);

	# Generate RNA-Seq QC Plot
	#my @METRICFILES;
	#for my $sampleName (sort keys %SAMPLES) {
	#	push @METRICSFILES, "analysis/$sampleName/$sampleName.rnaseqMetrics.txt";
	#}
	#$_ = join(" ", @METRICFILES);
	#`Rscript $dirname/multiRnaSeqCoverage.R $_ analysis/rnaseq.txt analysis/rnaseq.png`;
}


sub main {
	my $subsample = 0;
	my $species = 'hg19ERCC';
	my $tmpdir = '/scratch';
	my $file = '';

	if (scalar(@ARGV) == 4) {
		$subsample = $ARGV[0];
		$species = $ARGV[1];
		$tmpdir = $ARGV[2];
		$file = $ARGV[3];
	} else {
		print STDERR "Usage: $0 <subsample> <species> <tmpdir> <samples.txt>\n";
		exit(-1);
	}

	readSamples($file);
	runSamples($subsample, $species, $tmpdir);
	waitForSamples();
	collectAndPrintMetrics();
}

main();
