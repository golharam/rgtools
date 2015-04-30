#!/apps/sys/bin/perl
# This script uses the output from rgtools::CalculateTargetRegionCoverage
use strict;
use warnings;

my %SAMPLES;
my @ORDERED_SAMPLES;
my @ORDERED_TARGETS;
my %TARGETS;

if (scalar(@ARGV) != 1) {
	print STDERR "Usage: $0 <samples.txt>\n";
	exit 0;
}

# 1.  Get list of samples
open(IN, "<$ARGV[0]") || die "Unable to open $ARGV[0]";
while (<IN>) {
	chomp;
	my @fields = split(/\t/);
	my $s = $fields[0];
	my $per_target_file = "analysis/coverage/$s.targetRegionCoverage.txt";
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

	my $lineIndex = 24;
	if ($lines[$lineIndex] !~ m/^Target/) {
		die "File ".$SAMPLES{$s}{'file'}." different than expected\n";
	}

	# if first time, collect a list of targets
	if ($have_targets == 0) {
		$lineIndex = 25;
	        while ($lineIndex < scalar(@lines)) {
        	        my @fields = split(/\t/, $lines[$lineIndex]);
                	# This field may not always be defined
			#my $target = $fields[0];
                	my $target = $fields[1]."_".$fields[2]."_".$fields[3];
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
	
	$lineIndex = 25;
	my $targetIndex = 0;
	while ($lineIndex < scalar(@lines)) {
		my @fields = split(/\t/, $lines[$lineIndex]);
		#my $target = $fields[0];
		my $target = $fields[1]."_".$fields[2]."_".$fields[3];

		# make sure target matches our list
		if ($target ne $ORDERED_TARGETS[$targetIndex]) {
			die "Targets don't match";
		}

		$SAMPLES{$s}{$target}{'Total Reads'} = $fields[6];
		$SAMPLES{$s}{$target}{'Mean Coverage'} = $fields[7];
		$SAMPLES{$s}{$target}{'AvgPerBaseCoverage'} = $fields[8];
		$SAMPLES{$s}{$target}{'Base10X'} = $fields[9];
		$SAMPLES{$s}{$target}{'PctBases10X'} = $fields[10];
		$SAMPLES{$s}{$target}{'Base20X'} = $fields[11];
		$SAMPLES{$s}{$target}{'PctBases20X'} = $fields[12];
		$SAMPLES{$s}{$target}{'Base50X'} = $fields[13];
		$SAMPLES{$s}{$target}{'PctBases50X'} = $fields[14];
		$SAMPLES{$s}{$target}{'Base100X'} = $fields[15];
		$SAMPLES{$s}{$target}{'PctBases100X'} = $fields[16];

		$lineIndex++;
		$targetIndex++;
	}
}

# 3.  Output header rows
print "\t\t\t\t\t";
for my $s (@ORDERED_SAMPLES) {
	print "\t";
	print join("\t", "|----------", "-------------", $s, '--------', '--------','--------','--------','--------','--------','--------','----------|');
}
print "\n";
print join("\t", "Target", "Chromosome", "Start", "End", "Strand", "Length");
for my $s (@ORDERED_SAMPLES) {
	print "\t";
	print join("\t", "Total Reads", "Mean Coverage", "AvgPerBaseCoverage", 'Bases@10X', 'PctBases@10X',
									       'Bases@20X', 'PctBases@20X',
									       'Bases@50X', 'PctBases@50X',
									       'Bases@100X', 'PctBases@100X');
}
print "\n";

# 4.  Print out target and sample info
for my $target (@ORDERED_TARGETS) {
	print join("\t", $target, $TARGETS{$target}{'Chromosome'}, $TARGETS{$target}{'Start'}, $TARGETS{$target}{'End'}, 
			 $TARGETS{$target}{'Strand'}, $TARGETS{$target}{'Length'});

	for my $s (@ORDERED_SAMPLES) {
		print "\t";
		print join("\t", $SAMPLES{$s}{$target}{'Total Reads'}, $SAMPLES{$s}{$target}{'Mean Coverage'}, 
			$SAMPLES{$s}{$target}{'AvgPerBaseCoverage'}, 
			$SAMPLES{$s}{$target}{'Base10X'}, $SAMPLES{$s}{$target}{'PctBases10X'},
			$SAMPLES{$s}{$target}{'Base20X'}, $SAMPLES{$s}{$target}{'PctBases20X'},
			$SAMPLES{$s}{$target}{'Base50X'}, $SAMPLES{$s}{$target}{'PctBases50X'},
			$SAMPLES{$s}{$target}{'Base100X'}, $SAMPLES{$s}{$target}{'PctBases100X'}
			);
	}
	print "\n";
}

