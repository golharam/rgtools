print "# Created by Ryan Golhar\n";
print "# file is used by my scripts to keep track of samples and file locations for pipeline\n";
print "# Basedir: ", `pwd`;
print "#Sample\tFastQ1\tFastQ2\n";
for my $dir (`ls -d FASTQ/Sample_*`) {
	chomp $dir;
	my $f1 = `ls $dir/*R1*.fastq.gz`;
	chomp $f1;
	my $f2 = `ls $dir/*R2*.fastq.gz`;
	chomp $f2;

	print "$dir\t$f1\t$f2\n";
}
