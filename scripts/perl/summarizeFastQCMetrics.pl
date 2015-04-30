#!/usr/bin/env perl
use strict;

my %RESULTS;
my @CATEGORIES = ('Basic Statistics', 'Per base sequence quality', 'Per tile sequence quality', 'Per sequence quality scores',
		  'Per base sequence content', 'Per sequence GC content', 'Per base N content', 'Sequence Length Distribution',
		  'Sequence Duplication Levels', 'Overrepresented sequences', 'Adapter Content', 'Kmer Content');

# Scan the directory and locate at *_fastqc.zip file.  Extract the summary.txt and log the information.
my $analysisDir = $ARGV[0];

for my $zipFile (`ls $analysisDir/*_fastqc.zip`) {
	chomp $zipFile;
	
	my $fastQCDir = $zipFile;
	$fastQCDir =~ s/.zip$//;

	my $sample = $fastQCDir;
	$sample =~ s/$analysisDir\///;
	$sample =~ s/_fastqc//;

	`unzip $zipFile -d $analysisDir` if (! -d $fastQCDir);

	open(SUMMARY, "<$fastQCDir/summary.txt") || die "Unable to open $fastQCDir/summary.txt\n";
	while (<SUMMARY>) {
		chomp;
		m/(\w+)\t(.+)\t/;
		my $state = $1;
		my $category = $2;
		$RESULTS{$sample}{$category} = $state;
	}
	close(SUMMARY);

	`rm -rf $fastQCDir`;
}

my @SAMPLES = sort keys %RESULTS;
print "\t", join("\t", @SAMPLES), "\n";

for my $category (@CATEGORIES) {
	print "$category";

	for my $sample (@SAMPLES) {
		print "\t", $RESULTS{$sample}{$category};
	}
	print "\n";
}
