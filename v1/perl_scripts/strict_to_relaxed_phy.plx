#!/usr/bi/perl
#Oct. 7, 2011 edit to also remove all seqs comprised of only undermined values, i.e. all dashes
#also only print if #taxa included is more than 0
#Oct. 6, 2011 by Terri Porter
#Script to convert strict phylip (for ProtTest) into relaxed phylip format for RAxML 7.0.4
#usage perl strict_to_relaxed_phy.plx

use strict;
use warnings;

#declare  var
my $i=0;
my $infile;
my $j=0;
my $line;
my $filename;
my $ntax;
my $nchar;
my $seq;
my $taxon;

#declare array
my @infiles;
my @in;

#declare hash
my %alignment;

@infiles = qx(ls | grep .phy);

while ($infiles[$i]) {
	$infile = $infiles[$i];
	chomp $infile;

	open (IN,"<",$infile) || die "Error cannot read $infile:$!\n";
	@in = <IN>;
	close IN;

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		$filename = $infile.".relaxed";
		open (OUT,">>",$filename) || die "Error cannot write to $filename: $!\n";
		
		if ($line =~ /^\d+/) {
			$line =~ /^(\d+)\s+(\d+)/;
			$ntax = $1;
			$nchar = $2;
		}
		elsif ($line =~ /^[A-Z]/) {
			$line =~ tr/_/ /;
			$line =~ /(\w{2})\s+(\S+)/;
			$taxon = $1;
			$seq = $2;
			if ($seq !~ /[A-Z]/) {
				$ntax--;
			}
			else {
				$alignment{$taxon} = $seq;
			}
		}
		$j++;
	}

	if ($ntax > 0) {
		print OUT "$ntax\t$nchar\n";
		while (my($key,$value) = each (%alignment)) {
			print OUT "$key\t$value\n";
		}
	}
	elsif ($ntax==0) {
		unlink($filename);
	}
	$j=0;
	close OUT;
	$i++;
	@in=();
	%alignment=();
	$seq=();
	$taxon=();
	$nchar=();
	$ntax=();
	$filename=();
}
$i=0;
