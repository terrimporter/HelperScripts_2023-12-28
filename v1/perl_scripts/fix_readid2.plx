#!/usr/bin/perl
#April 5, 2013 by Terri Porter
#Script to fix truncated readids with full length readids in readid_taxonid.parsed2 lineage table (from taxonomy_crawl_map2.plx)
#The full length reference is 12F.fasta.sorted and the partial is 12F.nosingletons.fasta
#usage perl fix_readid.plx 12F.fasta.sorted readid_taxonid.parsed2

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $readid;
my $filename;
my $flag=0;
my $readidFrag;
my $old;
my $new;
my $fixreadid;
my $readidLine;
my $scalar;
my $readidPart;
my $newline;

#declare array
my @ref;
my @fix;
my @readidLine;
my @line;

#declare hash
my %readid; #indexed by full readid

open (REF, "<", $ARGV[0]) || die "Error cannot open reference file.fasta.sorted: $!\n";
@ref = <REF>;
close REF;

open (FIX, "<", $ARGV[1]) || die "Error cannot open file to be fixed file.nosingletons.fasta: $!\n";
@fix = <FIX>;
close FIX;

#hash full readids
while ($ref[$i]) {
	$line = $ref[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$readid = $1;
		$readidFrag = substr($readid,0,14);
		
		if (exists $readid{$readidFrag}) {
			$old = $readid{$readidFrag};
			$new = $old."|".$readid;
			$readid{$readidFrag} = $new;
		}
		else {
			$readid{$readidFrag} = $readid;
		}
	}
	$i++;
	$line=();
	$readid=();
	$readidFrag=();
	$old=();
	$new=();
}
$i=0;

$filename = $ARGV[1].".readidFixed";
open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";

#try to fix truncated readids
while ($fix[$i]) {
	$line = $fix[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$readidPart = $line[0];

	if ($readidPart =~ /^\d+/) {

		if (exists $readid{$readidPart}) {
			$readidLine = $readid{$readidPart};
			@readidLine = split(/\|/, $readidLine);
			$scalar = scalar(@readidLine);

			if ($scalar == 1) {
				$line[0] = $readidLine;
				$newline = join("\t", @line);
				print OUT "$newline\n";
			}
			else {
				print "PROBLEM: $readidPart maps to $readidLine\n";
				$flag=1;
			}
		}
		else {
			print "PROBLEM: can't find full readid for $readidPart\n";
			$flag=1;
		}
	}
	$flag=0;
	$i++;
	$line=();
	$fixreadid=();
	$readidLine=();
	@readidLine=();
	$scalar=();
	@line=();
	$readidPart=();
	$newline=();
}
$i=0;
close OUT;
