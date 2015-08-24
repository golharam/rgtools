#!/usr/bin/perl

# 1.  Make sure each gene exists only 1 one chr
my %GENE_TO_CHR;
open(GTF, "<$ARGV[0]") || die "Unable to open $ARGV[0]\n";
while (<GTF>) {
	chomp;
	my @fields = split(/\t/,$_);
	my $chr = $fields[0];
	my $attr = $fields[8];
	my $gene_id;
	if ($attr =~ m/gene_id \"(.+)\";/) {
		$gene_id = $1;
		if (defined($GENE_TO_CHR{$gene_id})) {
			if ($GENE_TO_CHR{$gene_id} ne $chr) {
				print STDERR "$gene_id exists on multiple chr: $GENE_TO_CHR{$gene_id}, $chr\n";
			}
		} else {
			$GENE_TO_CHR{$gene_id} = $chr;
		}
	} else {
		print STDERR "Unable to find gene_id for $attr\n";
	}
}
close(GTF);
