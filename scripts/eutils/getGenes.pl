#!/apps/sys/perl/bin/perl
use strict;
use LWP::Simple;

my $db = "gene";
my $query = "human[ORGN] AND alive[property]";

my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

my $output = get($url);

#parse WebEnv, QueryKey and Count
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

# retrieve data in batches of 20
my $retmax = 20;
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
	my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
	$efetch_url .= "&query_key=$key&retstart=$retstart";
	$efetch_url .= "&retmax=$retmax&retmode=xml";
	my $efetch_out = get($efetch_url);
	print "$efetch_out"
}
