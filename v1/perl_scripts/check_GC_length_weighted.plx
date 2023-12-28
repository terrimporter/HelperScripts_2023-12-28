#!/usr/bin/perl
#edited to weight by OTU size
#May 15, 2013 by Terri Porter
#Script to calculate %GC and length for sequences, do this separately for ITS1 and ITS2
#plot these values in excel/R to ensure correlation before pooling ITS1 and ITS2 for each site x temp for ITS1F
#usage perl check_GC_length.plx

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $flag=0;
my $header;
my $original;
my $new;
my $seq;
my $nextline;
my $j;
my $length;
my $base;
my $count=0;
my $minlength = 100; ### reset minimum length cutoff here ###
my $percentGC;
my $result;
my $abundance;

#declare array
my @ITS1;
my @ITS2;
my @seq;
my @header;

#declare hash
my %seq; #indexed by header
my %abundance; #indexed by header

open (ITS1, "<", "ITS1_99otus_minsize2.fa") || die "Error cannot open ITS1 file: $!\n";
@ITS1 = <ITS1>;
close ITS1;

open (ITS2, "<", "ITS2_99otus_minsize2.fa") || die "Error cannot open ITS2 file: $!\n";
@ITS2 = <ITS2>;
close ITS2;

#process ITS1
while ($ITS1[$i]) {
	$line = $ITS1[$i];
	chomp $line;

	if ($flag==0 && $line =~ /^>/) {
		$header = $line;
		@header = split(/;/,$header);
		$abundance = $header[1];
		$abundance =~ s/size=//g;
		$abundance{$header} = $abundance;
		$flag=1;
		$i++;
		next;
	}
	
	if ($flag==1 && $line !~ /^>/) { #accomodates seqs in fasta file divided by new lines
		$seq = $line;

		if (exists $seq{$header}) {
			$original = $seq{$header};
			$new = $original.$line;
			$seq{$header} = $new;
		}
		else {	
			$seq{$header} = $seq;
		}

		$j = $i+1;
		$nextline = $ITS1[$j];
		chomp $nextline;

		if ($nextline =~ /^>/) {
			$flag=0;
		}
		$i++;
		next;
	}
}
$i=0;

#parse through seqs in hash to determine length and GC%
#open (LENGTH, ">>", "ITS1.length") || die "Error cannot open ITS1 length outfile: $!\n";

#open (OUT, ">>", "ITS1_otus_min100bp.fa") || die "Error cannot open ITS1 otu outfile: $!\n";

open (GC, ">>", "ITS1_weighted.gc") || die "Error cannot open ITS1 gc outfile: $!\n";

while ( ($header,$seq) = each(%seq) ) {
	@seq = split('',$seq);
	$length = scalar(@seq);
#	print LENGTH "$length\n";

	if ($length >= $minlength) {
#		print LENGTH "$length\n";
#		print OUT ">$header\n$seq\n";
			
		while($seq[$i]) {
			$base = $seq[$i];

			if ($base =~ /(G|C)/) {
				$count++;
			}
			$i++;
			$base=();
		}
		$i=0;
		$percentGC = ($count/$length)*100;
		$abundance = $abundance{$header};
		$j=0;
		while ($j < $abundance) {
			print GC "$percentGC\n";
			$j++;
		}
		$count=0;
	}
	@seq=();
	$length=();
	$percentGC=();
	$abundance=();
}
#close LENGTH;
close OUT;
close GC;

$i=0;
$line=();
$header=();
$flag=0;
$seq=();
%seq=();
$original=();
$new=();
$j=();
$nextline=();
@seq=();
$length=();
$base=();
$count=0;
$percentGC=();
$abundance=();
%abundance=();

#process ITS2	
while ($ITS2[$i]) {
	$line = $ITS2[$i];
	chomp $line;

	if ($flag==0 && $line =~ /^>/) {
		$header = $line;
		@header = split(/;/,$header);
		$abundance = $header[1];
		$abundance =~ s/size=//g;
		$abundance{$header} = $abundance;
		$flag=1;
		$i++;
		next;
	}
	
	if ($flag==1 && $line !~ /^>/) { #accomodates seqs in fasta file divided by new lines
		$seq = $line;

		if (exists $seq{$header}) {
			$original = $seq{$header};
			$new = $original.$line;
			$seq{$header} = $new;
		}
		else {	
			$seq{$header} = $seq;
		}

		$j = $i+1;
		$nextline = $ITS2[$j];
		chomp $nextline;

		if ($nextline =~ /^>/) {
			$flag=0;
		}
		$i++;
		next;
	}
}
$i=0;

#parse through seqs in hash to determine length and GC%
#open (OUT, ">>", "ITS2_otus_min100bp.fa") || die "Error cannot open ITS2 otu outfile: $!\n";

open (GC, ">>", "ITS2_weighted.gc") || die "Error cannot open ITS2 gc outfile: $!\n";

while ( ($header,$seq) = each(%seq) ) {
	@seq = split('',$seq);
	$length = scalar(@seq);

	if ($length >= $minlength) {
#		print OUT "$header\n$seq\n";
			
		while($seq[$i]) {
			$base = $seq[$i];

			if ($base =~ /(G|C)/) {
				$count++;
			}
			$i++;
			$base=();
		}
		$i=0;
		$percentGC = ($count/$length)*100;
		$abundance = $abundance{$header};
		$j=0;
		while ($j < $abundance) {
			print GC "$percentGC\n";
			$j++;
		}
		$count=0;
	}
	@seq=();
	$length=();
	$percentGC=();
}
#close LENGTH;
close OUT;
close GC;

#$result = qx(grep ">" ITS1_otus_min100bp.fa | wc -l);
#chomp $result;
#print "$result\n";

#$result = qx(grep ">" ITS2_otus_min100bp.fa | wc -l);
#chomp $result;
#print "$result\n";

$result = qx(cat ITS1_weighted.gc ITS2_weighted.gc > ITS1ITS2_weighted.gc); #use this one for R?
