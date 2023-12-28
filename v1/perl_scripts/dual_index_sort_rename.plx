#!/usr/bin/perl
#Oct. 11, 2013 by Terri Porter
#Script to sort and rename reads according to dual indexes for Shadi's Barcoding Illumina Malaise dataset
#USAGE perl dual_index_sort_rename.plx fwd_tag.map rev_tag.map primers.txt derep.fasta

use warnings;
use strict;

#declare var
my $line;
my $i=0;
my $column; # 1-12
my $tag;
my $row; # A-H
my $name; # 'FWD' or 'REV_RC' only
my $primer;
my $header; 
my $k;
my $nextline;
my $potential_tag;
my $col_header;
my $filename;
my $prefix;
my $outfile;

#declare array
my @fwd;
my @rev;
my @primers;
my @derep;
my @line;
my @nextline;
my @temp;
my @filename;

#declare hash
my %fwd; #indexed by 5bp tag, value equals 96 well plate column number 1-12
my %rev; #indexed by 5bp tag, value equals 96 well plate row number A-H
my %primers; #indexed by 'FWD' or 'REV_RC', value equals primer sequence

### TAG + FWDprimer ###
### REV_RCprimer + TAG ###

open (FWD, "<", $ARGV[0]) || die "Error cannot open fwd_tag.map: $!\n";
@fwd = <FWD>;
close FWD;

open (REV, "<", $ARGV[1]) || die "Error cannot open rev_tag.map: $!\n";
@rev = <REV>;
close REV;

open (PRIMERS, "<", $ARGV[2]) || die "Error cannot open primers.txt: $!\n";
@primers = <PRIMERS>;
close PRIMERS;

open (DEREP, "<", $ARGV[3]) || die "Error cannot open derep.fasta: $!\n";
@derep = <DEREP>;
close DEREP;

$filename = $ARGV[3];
@filename = split(/\./,$filename);
$prefix = $filename[0];

#hash fwd tags
while ($fwd[$i]) {
	$line = $fwd[$i];
	chomp $line;
	
	if (length($line)>0) {
		@line = split(/\t/,$line);
		$column = $line[0];
		$tag = $line[1];
		$fwd{$tag} = $column;
	}

	$i++;
	$line=();
	@line=();
	$column=();
	$tag=();
}
$i=0;

#hash rev tags
while ($rev[$i]) {
	$line = $rev[$i];
	chomp $line;

	if (length($line)>0) {
		@line = split(/\t/,$line);
		$row = $line[0];
		$tag = $line[1];
		$rev{$tag} = $row;
	}

	$i++;
	$line=();
	@line=();
	$row=();
	$tag=();
}
$i=0;

#hash primers
while($primers[$i]) {
	$line = $primers[$i];
	chomp $line;

	if (length($line)>0) {
		@line = split(/\t/,$line);
		$name = $line[0];
		$primer = $line[1];
		$primers{$name} = $primer;
	}
	$i++;
	$line=();
	@line=();
	$name=();
	$primer=();
}
$i=0;

#check for degenerate bases in primers
#if present, save primer sequence as regex
while (($name,$primer) = each (%primers)) {
	if ($primer =~ /(B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z)/i) { #i.e. non-nucleotide bases I, N, Y, R, etc.
		reformat_degenerate_for_grep();
		$primers{$name} = $primer;
	}
	else {
		next;
	}
}

#search fasta file for TAG+FWDprimer, renaming as you go
open (OUT, ">>", "tempfile") || die "Error cannot open outfile: $!\n";

while ($derep[$i]) {
	$line = $derep[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(.+)/g;
		$header = $1;
		$k=$i+1;
		$nextline = $derep[$k];
		chomp $nextline;

		if (exists $primers{'FWD'}) {
			$primer = $primers{'FWD'};
			$nextline =~ /^(\w+)$primer/;
			$potential_tag = $1;

			if (exists $fwd{$potential_tag}) {
				$column = $fwd{$potential_tag};
				print OUT ">$column $header\n$nextline\n";
			}
			else {
				print "Error, cannot find potential tag $potential_tag in fwd hash\n";
			}
		}
		else {
			print "Error, could not find fwd primer in hash\n";
		}
		$i+=2;
	}
	else {
		$i++;
	}
	$line=();
	$header=();
	$k=();
	$nextline=();
	$primer=();
	$potential_tag=();
	$column=();
}
$i=0;

close OUT;

#search fasta file for REV_RCprimer+TAG, renaming as you go

open (TEMP, "<", "tempfile") || die "Error cannot open tempfile: $!\n";
@temp=<TEMP>;
close TEMP;

$outfile = $prefix.".derep.renamed.fasta";
open (OUT2, ">>", $outfile) || die "Error cannot open outfile: $!\n";

while ($temp[$i]) {
	$line = $temp[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(.+)/g;
		$col_header = $1;
		$col_header =~ /^(\d+)\s{1}(.+)/g;
		$column = $1;
		$header = $2;
		$k=$i+1;
		$nextline = $temp[$k];
		chomp $nextline;

		if (exists $primers{'REV_RC'}) {
			$primer = $primers{'REV_RC'};
			$nextline =~ /$primer(\w+)$/;
			$potential_tag = $1;

			if (exists $rev{$potential_tag}) {
				$row = $rev{$potential_tag};
				print OUT2 ">$column$row $header\n$nextline\n";
			}
			else {
				print "Error, cannot find potential tag $potential_tag in rev hash\n";
			}
		}
		else {
			print "Error, could not find rev primer in hash\n";
		}
		$i+=2;
	}
	else {
		$i++;
	}
	$line=();
	$header=();
	$k=();
	$nextline=();
	$primer=();
	$potential_tag=();
	$row=();
	$column=();
	$col_header=();
}
$i=0;

close OUT2;
unlink("tempfile");

####################

sub reformat_degenerate_for_grep {

my @bases;
my $base;
my $j=0;
my $newbase;

@bases = split(//,$primer);

while ($bases[$j]) {
	$base = $bases[$j];

	if ($base eq "U") {
		$newbase = "T";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "R") {
		$newbase = "[A|G]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "Y") {
		$newbase = "[C|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "S") {
		$newbase = "[G|C]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "W") {
		$newbase = "[A|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "K") {
		$newbase = "[G|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "M") {
		$newbase = "[A|C]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "B") {
		$newbase = "[C|G|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "D") {
		$newbase = "[A|G|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "H") {
		$newbase = "[A|C|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "V") {
		$newbase = "[A|C|G]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "N") { #anybase
		$newbase = "[A|C|G|T]";
		$bases[$j] = $newbase;
	}
	elsif ($base eq "I") { #inosine pairs with [ATC], Shadi used as if it's a "N"
		$newbase = "[A|C|G|T]";
		$bases[$j] = $newbase;
	}
	$j++;
	$base=();
	$newbase=();
}
$j=0;

$primer = join('', @bases)

}
