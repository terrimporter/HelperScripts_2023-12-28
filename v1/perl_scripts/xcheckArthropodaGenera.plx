#!/usr/bin/perl
#Terri Porter, Sept. 9, 2016
#Cross-check target Arthropoda genera (or whatever taxon) against Arthropoda v2 RDP classifier training set (testNBC.fasta and testNBC.taxonomy)
#USAGE perl xchecktargetgenera.plx testNBC.taxonomy testNBC.fasta targetArthropodaGenera.sortedunique

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $genus;
my $acc_kingdom;
my $acc;
my $list;
my $updatedlist;
my $target;
my $value;
my $flag=0;
my $count;

#declare array
my @taxonomy;
my @fasta;
my @targets;
my @line;
my @list;

#declare hash
my %genera_tax; #keys are genera, values are "1"
my %genera_fasta; #keys are genera, values are pipe delimited GenBank accessions

open (TAX, "<", $ARGV[0]) || die "Cannot open first infile: $!\n";
@taxonomy = <TAX>;
close TAX;

open (FASTA, "<", $ARGV[1]) || die "Cannot open second infile: $!\n";
@fasta = <FASTA>;
close FASTA;

open (TARGET, "<", $ARGV[2]) || die "Cannot open third infile: $!\n";
@targets = <TARGET>;
close TARGET;

#parse testNBC.taxonomy into a hash with genus keys
while ($taxonomy[$i]) {
	$line = $taxonomy[$i];
	chomp $line;

	if ($line =~ /genus$/) {
		@line = split(/\*/, $line);
		$genus = $line[1]; #genus[1]
		$genera_tax{$genus}="1";
	}
	$i++;
	$line=();
	@line=();
	$genus=();
}
$i=0;

#parse testNBC.fasta into a hash with genus keys and pipe delimited list of GenBank accessions
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		@line = split(/;/, $line);
		$acc_kingdom = $line[0]; #>GenbankAccession Metazoa
		$acc_kingdom =~ s/^>//g;
		$acc_kingdom =~ s/ Metazoa//g;
		$acc = $acc_kingdom;
		$genus = $line[5]; #genus

		if ($genera_fasta{$genus}) {
			$list = $genera_fasta{$genus};
			$updatedlist = $list."|".$acc;
			$genera_fasta{$genus} = $updatedlist;
		}
		else {
			$genera_fasta{$genus} = $acc;
		}
	}
	$i++;
	$line=();
	@line=();
	$acc_kingdom=();
	$acc=();
	$genus=();
	$list=();
	$updatedlist=();

}
$i=0;

open (OUT, ">>", "target_arthropoda_in_reference.txt") || die "Cannot open outfile:$!\n";
print OUT "TargetArthropodaGenera\tPresentInArthropodaTaxonomy\tNumInArthropodaFasta\n";

#for each target genus, check both hashes, record presence in testNBC.taxonomy and number of accessions in testNBC.fasta
while ($targets[$i]) {
	$target = $targets[$i];
	chomp $target;

	print OUT "$target\t";

	if (exists $genera_tax{$target}) {
		$value = $genera_tax{$target};
		print OUT "$value\t";
		$flag=1;
	}
	else {
		print OUT "0\t0\n";
	}

	if ($flag==1) {
		if (exists $genera_fasta{$target}) {
			$list = $genera_fasta{$target};

			if ($list =~ /\|/) {
				@list = split(/\|/, $list);	
				$count = scalar(@list);

				print OUT "$count\n";
			}
			else {
				print OUT "1\n";
			}
		}
		else {
			print "Error, target genus $target found in testNBC.taxonomy but not in testNBC.fasta\n";
		}
	}

	$i++;
	$target=();
	$value=();
	$flag=0;
	$list=();
	@list=();
	$count=();
}
$i=0;
close OUT;
