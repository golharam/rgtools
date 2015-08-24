#!/usr/bin/env perl
#use strict;
# This script reports nonsynonymous variants in Xpress format

use Getopt::Long;
use Vcf;

my $vcfFile;
GetOptions("vcf|v=s" => \$vcfFile) || die "Error in command link arguments\n";

my $vcf = Vcf->new(file=>$vcfFile);
$vcf->parse_header();
my (@samples) = $vcf->get_samples();
my @variants;
my @DATA;

# Print the header ie variants we are interested in reporting
while (my $x=$vcf->next_data_hash()) {
	if ($x->{'INFO'}{'SNPEFF_EFFECT'} =~ m/NON_SYNONYMOUS_CODING|START_GAINED/) {
		my $genename = $x->{'INFO'}{'SNPEFF_GENE_NAME'};
		if ($x->{'INFO'}{'SNPEFF_EFFECT'} =~ m/NON_SYNONYMOUS_CODING/) {
			$genename .= '('.$x->{'INFO'}{'SNPEFF_AMINO_ACID_CHANGE'}.')';
		} else {
			$genename .= '('.$x->{'INFO'}{'SNPEFF_EFFECT'}.')';
		}
		print "\t$genename";

		push @variants, $x;
	}
}
print "\n";

# Print the sample info per variant
for my $sample (@samples) {
	print "$sample";

	for my $variant (@variants) {
		print "\t";
		print $variant->{'gtypes'}{$sample}{'GT'} if ($variant->{'gtypes'}{$sample}{'GT'} ne '0/0');
	}
	print "\n";
}
