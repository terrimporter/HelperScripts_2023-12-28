#!/usr/bin/perl
#June 12, 2012 by Terri Porter
#Script to filter out contam sequences from seed.fasta.list (from make_mothur_list_file.plx)
#usage perl filter_contam_from_mothur_list.plx seed.fasta.list 12F.contam.list

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $j=0;
my $part;
my $flag=0;
my $contam;
my $original;
my $partcontam;
my $k=0;

#declare array
my @in;
my @line;
my @contam;
my @part;

#declare hash
my %contam;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (CONTAM, "<", $ARGV[1]) || die "Error cannot open infile2: $!\n";
@contam = <CONTAM>;
close CONTAM;

#add comtams ids into a form for regex
while ($contam[$i]) {
	$line = $contam[$i];
	chomp $line;

	if ($i==0) {
		$contam = "(".$line;
	}
	else {
		$original = $contam;
		$contam = $original."|".$line;
	}
	$i++;
	$line=();
	$original=();
}
$i=0;
$contam = $contam.")";
print "contam: $contam\n";

open (OUT, ">>", "seed.fasta.list2") || die "Error cannot open outfile: $!\n";

open (OUT2,">>", "contam.full") || die "Error cannot open contam full outfile: $!\n";

#parse through the mothur list
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);

	while ($line[$j]) {
		$part = $line[$j];
		
		if ($part =~ /usearch/) {
			print OUT "$part\t";
			$flag=1;
		}
		elsif ($flag==1) {
			print OUT "$part\t";
			$flag=0;
		}
		elsif ($part =~ /$contam/) {
			@part = split(',', $part);
			
			while ($part[$k]) {
				$partcontam = $part[$k];
				print OUT2 "$partcontam\n";
				$k++;
			}
			$k=0;

			$j++;
			next;
		}
		else {
			print OUT $part."\t";
		}
		$j++;
		$part=();
		$flag=0;
	}
	$j=0;
	$i++;
	$line=();
	@line=();
}
$i=0;
close OUT;
close OUT2;
