#!/usr/bin/perl
use List::MoreUtils qw(uniq);
use strict;

my %ACCtoSYM;

# Load gene2refseq for Accession to Symbol
open(IN, "<gene2refseq") || die "Unable to open gene2refseq";
while(<IN>) {
	next if ($_ !~ m/^9606/);
	
	chomp;
	my @fields = split(/\t/);
	$fields[3] =~ m/(\w\w\_\d+)/;
	my $accession = $1;
	my $genesymbol = $fields[15];
	$ACCtoSYM{$accession} = $genesymbol;
}
close(IN);
print "Loaded ", scalar (keys %ACCtoSYM), " accessions mapped to ", scalar uniq((values %ACCtoSYM)), " gene symbols.\n";

open(GTF, "<$ARGV[0]") || die "Unable to open GTF File: $ARGV[0]\n";
while (<GTF>) {
	chomp;
	my @fields = split(/\t/);
	my $attr = $fields[8];
 	
 	my $gene_id;
        if ($attr =~ m/gene_id \"(.+)\";/) {
                $gene_id = $1;
        } else {
                print STDERR "Unable to find gene_id in $attr\n";
        }

 	my $transcript_id;
        if ($attr =~ m/transcript_id \"(.+)\";/) {
                $transcript_id = $1;
        } else {
                print STDERR "Unable to find transcript_id in $attr\n";
        }

	if (!defined($ACCtoSYM{$transcript_id})) {
		print STDERR "No map for $transcript_id\n"
	}
	
}
close(GTF);
