#!/usr/bin/perl
#Aug. 29, 2012 by Terri Porter
#Script to parse a directory of fasta files from estOrtholgos for each ortho, look for yeast seqs, and grab orthoID
#usage perl grab_yeast_GOid.plx

use strict;
use warnings;

#declare var
my $dir;
my $i=0;
my $file;
my $pathToFile;
my $j=0;
my $line;
my $id;

#declare array
my @files;
my @in;
my @ids;

print "Enter path containing ortho fasta files from estOrthologs.plx including final slash:\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error cannot open directory: $!\n";
@files = readdir(DIR);
closedir (DIR);

while ($files[$i]) {
	$file = $files[$i];
	chomp $file;

	$pathToFile = $dir.$file;

	open(IN,"<",$pathToFile) || die "Error cannot open $pathToFile: $!\n";
	@in = <IN>;
	close IN;

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($line =~ /^>SC/) {
			$line =~ /^>SC(\d+),/;
			$id = $1;
			push(@ids, $id);
		}
		$j++;

		$line=();
		$id=();
	}
	$i++;

	$j=0;
	$file=();
	$pathToFile=();
	@in=();
}
$i=0;

open (OUT,">>", "yeast.ids") || die "Error cannot open outfile: $!\n";

while ($ids[$i]) {
	$id = $ids[$i];
	chomp $id;

	print OUT "$id\n";
	$i++;
}
$i=0;

close OUT;
