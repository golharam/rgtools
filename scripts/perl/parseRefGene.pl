#!/apps/sys/perl/bin/perl
use strict;

open(IN, "<$ARGV[0]") || die "Unable to open $ARGV[0]\n";
$_ = <IN>;
my @header = split(/\t/, $_);

my %GENES;
my %ACC2GENE;

while (<IN>) {
	chomp;
	my @fields = split(/\t/, $_);
	if ($#header != $#fields) {
		print STDERR "Found incorrect number of fields\n";
	} else {
		my $acc = $fields[1];
		my $gene = $fields[12];
		$GENES{$gene}{$acc} += 1;
		$ACC2GENE{$acc} = $gene;
	}
}
close(IN);

open(IN, "<$ARGV[1]") || die "Unable to open GTF File, $ARGV[1]\n";
while (<IN>) {
	chomp;
	my @fields = split(/\t/, $_);
	$fields[8] =~ /gene_id \"(NM_\d+)\"/;
	my $gene = $ACC2GENE{$1};
	$fields[8] =~ s/gene_id \"(NM_\d+)\"/gene_id \"$gene\"/;

	print join("\t", @fields), "\n";
}
close(IN);
#for my $gene (sort keys %GENES) {
#	print "$gene";

#	for my $acc (keys $GENES{$gene}) {
#		print "\t$acc";
#	}
#	print "\n";
#}
