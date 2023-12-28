#!/usr/bin/perl
#Oct. 16, 2013, tried to use agrep with regex but ignored indels and substitutions, so avoid looking for degenerate primers and just search for tags for now
#Oct. 11, 2013 by Terri Porter
#Script to sort and rename reads according to dual indexes for Shadi's Barcoding Illumina Malaise dataset
#USAGE perl dual_index_sort_rename.plx fwd_tag.map rev_tag.map derep.fasta

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
my $header1;
my $header2;

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

open (DEREP, "<", $ARGV[2]) || die "Error cannot open derep.fasta: $!\n";
@derep = <DEREP>;
close DEREP;

$filename = $ARGV[2];
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

#search fasta file for forward tag, renaming as you go
open (OUT, ">>", "tempfile") || die "Error cannot open outfile: $!\n";

$i=0;
while ($derep[$i]) {
	$line = $derep[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(.+)$/g;
		$header=$1;
		$k=$i+1;
		$nextline = $derep[$k];
		chomp $nextline;

		if ($nextline =~ /^(\w{5})/) {
			$nextline =~ /^(\w{5})/;
			$potential_tag = $1;
#			print "$potential_tag\n";

			if (exists $fwd{$potential_tag}) {
				$column = $fwd{$potential_tag};
				print OUT ">$column $header\n$nextline\n";
			}
			else {
#				print "Fasta file: Line $k+1 - cannot find potential tag $potential_tag in fwd hash\n";
			}
		}
		else {
#			print "Fasta file: Line $k+1 - no tag match\n";
		}
	}
	$i++;
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

		if ($col_header =~ /^(\d+)\s{1}(.+)/g) {
			$col_header =~ /^(\d+)\s{1}(.+)/g;
			$column = $1;
			$header = $2;
			$k=$i+1;
			$nextline = $temp[$k];
			chomp $nextline;

			if ($nextline =~ /(\w{5})$/) {
				$nextline =~ /(\w{5})$/;
				$potential_tag = $1;

				if (exists $rev{$potential_tag}) {
					$row = $rev{$potential_tag};
					print OUT2 ">$column$row $header\n$nextline\n";
				}
				else {
#					print "Tempfile $k - cannot find potential tag $potential_tag in rev hash\n";
				}
			}
			else {
#				print "Tempfile Line $k+1 - No good tag match\n";
			}		
		}
		else {
#			print "Tempfile Line $k+1 - cannot find column in header\n";
		}
	}
	$i++;
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


