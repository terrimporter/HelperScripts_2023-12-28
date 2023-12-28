#!/usr/bin/perl
#Oct. 11, 2013 by Terri Porter
#Script to hash the output from dual_index_sort_rename.plx for Shadi's Barcoding illumina malaise dataset
#usage perl hash_renamed_fastas.plx prefix.derep.renamed.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $plateloc; #column (1-12) - row (A-H) location on a 96 well plate
my $j;
my $nextline;
my $original;
my $new;
my $seqlist;
my $scalar;
my $count;
my $filename;
my $prefix;
my $outfile;

#declare array
my @in;
my @seqlist;
my @filename;

#declare hash
my %seq; #indexed by plateloc, value eq list of sequences


open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

$filename = $ARGV[0];
@filename = split(/\./,$filename);
$prefix = $filename[0];

#hash infile
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\d+\w+)\s{1}.+$/g;
		$plateloc = $1;
		$j = $i+1;
		$nextline = $in[$j];
		chomp $nextline;
#		print $plateloc."\n";
		if (exists $seq{$plateloc}) {
			$original = $seq{$plateloc};
			$new = $original."|".$nextline;
			$seq{$plateloc} = $new;
		}
		else {
			$seq{$plateloc} = $nextline;
		}
	}
	$i++;
	$line=();
	$plateloc=();
	$j=();
	$nextline=();
	$original=();
	$new=();
}
$i=0;

#count %seq hash keys to see how many individuals were successfully dual-tagged and recovered
$count = keys %seq;
print "$count individuals were successfully dual-tagged and recovered\n";

#print each individual plateloc to a file along with seq count
$outfile = $prefix.".recovered.txt";
open (OUT, ">>", $outfile) || "Error cannot open outfile:$!\n";

while (($plateloc,$seqlist) = each (%seq)) {
	@seqlist = split(/\|/,$seqlist);
	$scalar = scalar(@seqlist);

	print OUT "$plateloc\t$scalar\n";

	@seqlist=();
	$scalar=();
}
close OUT;
