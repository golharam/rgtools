#!/apps/sys/perl/bin/perl
use strict;

open(GTF, "<$ARGV[0]") || die "Unable to open $ARGV[0]\n";
while (<GTF>) {
	if (m/^#/) {
		print $_;
		next;
	}
	print 'chr'.$_;
}
close(GTF);
