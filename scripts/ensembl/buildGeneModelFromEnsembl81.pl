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

v2($GTF);

sub v2 {
	my ($GTF) = @_;
	
	# From http://crazyhottommy.blogspot.com/2013/05/find-exons-introns-and-intergenic.html
	
	`cat $GTF | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4,$5}' > exons.bed`;
	`sortBed -i exons.bed > exons.sorted.bed`;
	# In mouse, there is a problem with one line that needs to be edited manually
	`mergeBed -i exons.sorted.bed > exons.bed`;  
	`rm exons.sorted.bed`;

	`cat $GTF | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5}' > genes.bed`;
	`sortBed -i genes.bed > genes.sorted.bed`;
	`mv genes.sorted.bed gene.bed`;

	`subtractBed -a genes.bed -b exons.bed > introns.bed`;

	`complementBed -i genes.bed -g ../../$SPECIES.fa.fai > intergenic.bed`;

	`grep "gene_biotype \"rRNA\"" $GTF | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5}' > rRNA.bed`;
	`sortBed -i rRNA.bed > rRNA.sorted.bed`;
	`mv rRNA.sorted.bed rRNA.bed`;

}

sub v1 {
	my ($GTF) = @_;
	print STDERR "Creating GTF files...\n";
	open(GTF, "<$GTF") || die "Unable to open $GTF\n";
	open(GENESGTF, ">genes.gtf") || die "Unable to open genes.gtf\n";
	open(EXONSGTF, ">exons.gtf") || die "Unable to open exons.gtf\n";
	open(RRNAGTF, ">rRNA.gtf") || die "Unable to open rRNA.gtf\n";
	while(<GTF>) {
		next if (m/^#/);
		my @fields = split(/\t/, $_);
	
		if (($fields[8] =~ m/gene_biotype "rRNA"/) && ($fields[2] =~ m/gene/)) {
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
	
	`mv genes.sorted.gtf genes.gtf`;
	`mv exons.sorted.gtf exons.gtf`;
	`mv introns.sorted.gtf introns.gtf`;
	`mv rRNA.sorted.gtf rRNA.gtf`;
	
	print STDERR "Indexing GTF files\n";
	`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index genes.gtf`;
	`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index exons.gtf`;
	`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index introns.gtf`;
	`/apps/sys/galaxy/external_packages/igvtools-2.3.57/IGVTools/igvtools index rRNA.gtf`;
	
	print STDERR "Convert GTF to BED\n";
	`gtf2bed < genes.gtf > genes.bed`;
	`gtf2bed < exons.gtf > exons.bed`;
	`gtf2bed < introns.gtf > introns.bed`;
	`gtf2bed < rRNA.gtf > rRNA.bed`;
