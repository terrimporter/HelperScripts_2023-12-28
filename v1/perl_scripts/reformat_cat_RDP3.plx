#!/usr/bin/perl
#Dec. 5, 2016 edited to parse Gibson et al., 2015 BE and F230 clustered reads from Dryad
#Terri Porter, Sept. 7, 2016; Script to add sample name and amplicon name to each RDP taxonomic assignment, concatenate all into a single file, import into excel
#USAGE perl reformat_cat_RDP.plx

use strict;
use warnings;

#declare var
my $dir;
my $i=0; #dir file counter
my $j=0; #file line counter
my $line;
my $sample; 
my $amplicon;
my $outfile;
my $file;
my $pathtofile;

#declare array
my @files;
my @contents;

print "Enter directory containing RDP out files (including final / ):\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error opening dir $dir\n";
@files = readdir DIR;
closedir DIR;

while ($files[$i]) {
	$file = $files[$i];
	chomp $file; ###i.e. format RR01.BE.fasta.out

	if ($file =~/^\w/) {
		print "file: $file\n";#test
		$sample = substr($file,0,3); #i.e. RR01
#	$amplicon = substr($file,5,2); #i.e. BE
		$pathtofile = $dir.$file;

		open (IN, "<", $pathtofile) || die "Error cannot open file $file: $!\n";
		@contents = <IN>;
		close IN;

#		$outfile = $dir."concatenated.".$amplicon.".txt";
		$outfile = $dir."concatenated.txt";

		open (OUT, ">>", $outfile) || die "Error cannot open outfile $outfile: $!\n";

		while ($contents[$j]) {
			$line = $contents[$j];
			chomp $line;

			if ($line=~/;size=/) { ####size- or size= change pattern where appropriate
#				$line =~s/;size=/;$sample;$amplicon;/g;
				$line =~s/BenthicPAD-BE-/BE\t/g;
				$line =~s/_S\d{2}_L001_PAIRED_GLGQ_001_Seq\d+;/\t/g;
#				$line =~s/;size=/;$sample;/g;
#				$line =~s/\t\t//g;
				$line =~s/superkingdom\t//g;
				$line =~s/kingdom\t//g;
				$line =~s/phylum\t//g;
				$line =~s/class\t//g;
				$line =~s/order\t//g;
				$line =~s/family\t//g;
				$line =~s/genus\t//g;
				print OUT $line."\n";
			}
			
			$j++;
		}
		$j=0;
		close OUT;
	}
	$i++;

}
$i=0;
