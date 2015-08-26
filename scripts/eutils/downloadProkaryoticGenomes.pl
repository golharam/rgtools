#!/bin/env perl
# Taken from http://www.ncbi.nlm.nih.gov/books/NBK25498/pdf/Bookshelf_NBK25498.pdf
# Application 2: Converting accession numbers to data
use strict;
use warnings;
use LWP::Simple;

my @REFSEQIDS = `tail -n +3 ProkayyoticReferenceGenomes.txt  | cut -f 4`;

# Since some lines contains multiple ids seperated by commas, 
# it easier to build a CSV list then split them up.
chomp @REFSEQIDS;
my $csvlist = join(',', @REFSEQIDS);
@REFSEQIDS = split(',', $csvlist);
print STDERR "Searching for ", scalar(@REFSEQIDS), " ids.\n";

# append [accn] field to each accession
for (my $i = 0; $i < @REFSEQIDS; $i++) {
	$REFSEQIDS[$i] .= '[accn]';
}

# join the accessions with OR
my $query = join('+OR+', @REFSEQIDS);

# assemble the esearch URL
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";

# post the esearch URL
my $output = get($url);

# parse WebEnv and QueryKey
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

print STDERR "Found: $count records\n";
if ($count == 0) {
    exit(0);
}

# open output file for writing
open(OUT, ">tmp.prokaryotes.fa") || die "Can't open file!\n";

# retrieve data in batches of 500
my $retmax = 500;
for (my $ret = 0; $ret < $count; ) {
    my $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
    $efetch_url .= "&query_key=$key&retstart=$ret";
    $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
    print STDERR "Fetching $ret - $retmax\n";

    my $efetch_out = get($efetch_url);
    my $actual_sequences_returned = $efetch_out =~ s/>/\n>/g;  # count number of sequences returned
    $ret += $actual_sequences_returned;
    print OUT "$efetch_out";
    print STDERR "Fetched $ret\n";
}
close OUT;

rename("tmp.prokaryotes.fa", "prokayotes.fa");

