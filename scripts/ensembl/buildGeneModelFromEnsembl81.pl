#!/apps/sys/perl/bin/perl
use strict;
use warnings;

# List of species available at:
# http://useast.ensembl.org/info/data/ftp/index.html
# This script currently works for Human Release 81

if (scalar(@ARGV) != 1) {
	die "Usage: $0 <Ensembl GTF File>\n";
}
my $GTF = shift @ARGV;

print STDERR "Creating GTF files...\n";
open(GTF, "<$GTF") || die "Unable to open $GTF\n";
open(GENESGTF, ">genes.gtf") || die "Unable to open genes.gtf\n";
open(EXONSGTF, ">exons.gtf") || die "Unable to open exons.gtf\n";
open(RRNAGTF, ">rRNA.gtf") || die "Unable to open rRNA.gtf\n";
while(<GTF>) {
	next if (m/^#/);
	my @fields = split(/\t/, $_);
	if (($fields[1] =~ m/rRNA/) && ($fields[2] =~ m/gene/)) {
		print RRNAGTF $_;
	} elsif ($fields[2] =~ m/gene/) {
		print GENESGTF $_;
	} elsif ($fields[2] =~ m/exon/) {
		print EXONSGTF $_;
	}
}
close(RRNAGTF);
close(EXONSGTF);
close(GENESGTF);
close(GTF);

print STDERR "Sorting GTF files\n";
`sortBed -i genes.gtf > genes.sorted.gtf`;
`sortBed -i exons.gtf > exons.sorted.gtf`;
`sortBed -i rRNA.gtf > rRNA.sorted.gtf`;
`subtractBed -s -a genes.gtf -b exons.gtf > introns.gtf`;
`sortBed -i introns.gtf > introns.sorted.gtf`;

print STDERR "Indexing GTF files\n";
`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index genes.sorted.gtf`;
`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index exons.sorted.gtf`;
`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index introns.sorted.gtf`;
`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index rRNA.sorted.gtf`;

print STDERR "Convert GTF to BED\n";
`gtf2bed < genes.sorted.gtf > genes.bed`;
`gtf2bed < exons.sorted.gtf > exons.bed`;
`gtf2bed < introns.sorted.gtf > introns.bed`;
`gtf2bed < rRNA.sorted.gtf > rRNA.bed`;
