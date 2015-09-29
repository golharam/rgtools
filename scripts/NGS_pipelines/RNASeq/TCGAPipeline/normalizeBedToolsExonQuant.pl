#!/usr/bin/env perl
$bamfile = $ARGV[0];
$exonBed = $ARGV[1];

my $command = "samtools view $bamfile | head -1 | cut -f 10";
my $read = `$command` or die $1;
my $readLength = length($read) - 1; # Subtract 1 to discard newline

$i=0;
open IN, $exonBed or die "cant open $exonBed: $!\n";
@outputOrder = ();
while(<IN>){
	chomp;
	(undef, undef, undef, undef, $id) = split /\t/;
	$outputOrder[$i++] = $id;
}
close IN;


$tot = 0;
while(<STDIN>){
	(undef, undef, undef, undef, $id, $counts, undef, $length, undef) = split /\t/;
	$tot += $counts;
	$exonDat{$id} = $counts . "\t" . $length;
}

for($i = 0; $i<=$#outputOrder; $i++){
	$id = $outputOrder[$i];
	next unless defined($exonDat{$id});
	($counts, $length) = split /\t/, $exonDat{$id};
	$cov = ($counts*$readLength)/($length+1);
	$rpkm = 10**9 * $counts/(($length+1)*$tot);
	print $id . "\t" . $counts . "\t" . $cov . "\t" . $rpkm . "\n";
} 
