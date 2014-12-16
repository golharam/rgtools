#!/apps/sys/bin/perl
# This script uses the output from rgtools::CalculateTargetRegionCoverage
use strict;
use warnings;

my %SAMPLES;
my @ORDERED_SAMPLES;
my @ORDERED_TARGETS;
my %TARGETS;

# 1.  Get list of samples
open(IN, "<$ARGV[0]") || die "Unable to open $ARGV[0]";
while (<IN>) {
	chomp;
	my @fields = split(/\t/);
	my $s = $fields[0];
	my $per_target_file = `ls $s.*.coverage_per_target_metrics.txt`;
	chomp $per_target_file;
	if (! -e $per_target_file) {
		die "Unable to find $per_target_file";
	}
	$SAMPLES{$s}{'file'} = $per_target_file;
	push @ORDERED_SAMPLES, $s;
}
close(IN);

# 2.  Get a list of targets
my $have_targets = 0;
for my $s (@ORDERED_SAMPLES) {
	open(IN, "<".$SAMPLES{$s}{'file'}) || die "Unable to open ".$SAMPLES{$s}{'file'};
	my @lines = <IN>;
	chomp @lines;
	close(IN);

	my $lineIndex = 21;
	if ($lines[$lineIndex] !~ m/^Target/) {
		die "File ".$SAMPLES{$s}{'file'}." different than expected\n";
	}

	# if first time, collect a list of targets
	if ($have_targets == 0) {
		$lineIndex = 22;
	        while ($lineIndex < scalar(@lines)) {
        	        my @fields = split(/\t/, $lines[$lineIndex]);
                	my $target = $fields[0];
	                push @ORDERED_TARGETS, $target;
        	        $TARGETS{$target}{'Chromosome'} = $fields[1];
                	$TARGETS{$target}{'Start'} = $fields[2];
	                $TARGETS{$target}{'End'} = $fields[3];
        	        $TARGETS{$target}{'Strand'} = $fields[4];
                	$TARGETS{$target}{'Length'} = $fields[5];
			$lineIndex++;
		}
		$have_targets = 1;
	}
	
	$lineIndex = 22;
	my $targetIndex = 0;
	while ($lineIndex < scalar(@lines)) {
		my @fields = split(/\t/, $lines[$lineIndex]);
		my $target = $fields[0];
		# make sure target matches our list
		if ($target ne $ORDERED_TARGETS[$targetIndex]) {
			die "Targets don't match";
		}

		$SAMPLES{$s}{$target}{'Total Reads'} = $fields[6];
		$SAMPLES{$s}{$target}{'Mean Coverage'} = $fields[7];
		$SAMPLES{$s}{$target}{'AvgPerBaseCoverage'} = $fields[8];
		$SAMPLES{$s}{$target}{'Base10X'} = $fields[9];
		$SAMPLES{$s}{$target}{'PctBases10X'} = $fields[10];

		$lineIndex++;
		$targetIndex++;
	}
}

# 3.  Output header rows
open(OUT, ">pertargetmetrics.txt") || die "Unable to open pertargetmetrics.txt";
print OUT "\t\t\t\t\t";
for my $s (@ORDERED_SAMPLES) {
	print OUT "\t";
	print OUT join("\t", "|----------", "-------------", $s, '--------', '----------|');
}
print OUT "\n";
print OUT join("\t", "Target", "Chromosome", "Start", "End", "Strand", "Length");
for my $s (@ORDERED_SAMPLES) {
	print OUT "\t";
	print OUT join("\t", "Total Reads", "Mean Coverage", "AvgPerBaseCoverage", 'Bases@10X', 'PctBases@10X');
}
print OUT "\n";

# 4.  Print out target and sample info
for my $target (@ORDERED_TARGETS) {
	print OUT join("\t", $target, $TARGETS{$target}{'Chromosome'}, $TARGETS{$target}{'Start'}, $TARGETS{$target}{'End'}, 
			 $TARGETS{$target}{'Strand'}, $TARGETS{$target}{'Length'});

	for my $s (@ORDERED_SAMPLES) {
		print OUT "\t";
		print OUT join("\t", $SAMPLES{$s}{$target}{'Total Reads'}, $SAMPLES{$s}{$target}{'Mean Coverage'}, 
			$SAMPLES{$s}{$target}{'AvgPerBaseCoverage'}, $SAMPLES{$s}{$target}{'Base10X'},
			$SAMPLES{$s}{$target}{'PctBases10X'});
	}
	print OUT "\n";
}
close(OUT);

