#!/apps/sys/perl/bin/perl
use strict;
use LWP::Simple;
use XML::LibXML;

my $id = '7157';
$id = $ARGV[0] if (scalar(@ARGV) == 1);

print "Retreiving gene id $id\n";
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "efetch.fcgi?db=gene&id=$id&retmode=xml";

my $out = get($url);

#print $out;
print "Parsing...\n";
my $parser = XML::LibXML->new();
my $doc = $parser->load_xml(string => $out);

my $genename = ($doc->findnodes("//Gene-ref_locus"))->to_literal;
print $genename, "\n";

my %ACCESSION;
foreach my $geneset_accession ($doc->findnodes('//Gene-commentary_accession')) {
	my $acc = $geneset_accession->to_literal;
	if ($acc =~ m/^NM_/) {
		$ACCESSION{$acc} = 1;
	}
}

for my $acc (keys %ACCESSION) {
	print "$acc\n";
}
