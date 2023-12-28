#!/usr/bin/perl
#January 27, 2010 by Terri Porter
#Updated Aug.31, 2010 to accept infile as an argument instead
#Script to parse genbank records retrieved using the script get_genbankrecord.plx and retrieve organism lineage.
#Usage #perl parse_genbankrecord.plx
#Script will automatically look for gb.outfile from get_genbankrecord.plx and will produce parsed_gb.outfile which is a tab-delimited file
#new usage $perl parse_genbankrecord_using_gi.plx file.gb
#March22, 2011 tweak to parse gb files individually and grab only the entry matching a specific gi, ignorning the rest of the entry.
#new usage: perl parse_genbankrecord_using_gi2.plx

use strict;
use warnings;

#declare variables
my $line;
my $a=0; #flags start of new record
my $gi;
my $gb;
my $b=0; #flags accession field
my $lineage;
my $c=0; #flags orgname field
my $d=0; #flags lineage field
my $k; #counter
my $flag=0;
my $dir;
my $i=2;
my $filename;
my $file;
my $path_to_file;

#declare arrays
my @record;
my @files;
my @file;

print "Enter full path to directory containing .gb files (including trailing '/'):\n";
$dir = <STDIN>;
chomp $dir;
#$dir = "$dir";

opendir DH, $dir;
@files = readdir (DH);
print "@files\n";#test - ok
while ($files[$i]) {
	$file = $files[$i];
	$path_to_file = $dir.$file;
	@file = split(/\./,$file);
	$gi = $file[0];
	
	open (IN, '<', $path_to_file) || die ("Error cannot read .gb file: $!\n");
	@record = <IN>;
	chomp(@record);
	close IN;

	$filename = $gi.".parsed";

	open (OUT, '>>', $filename) || die ("Error cannot write .parsed: $!\n");

	foreach (@record) {
		$line = $_;
		chomp $line;
		
		if ($c==0) {

		if ($a==0) {
			if ($line =~ /$gi/) {
				$a = 1;
				$k++; #entry counter
				next;
			}
		}
		elsif ($a==1) {
			if ($flag ==0 ) {
				$line =~ /descr\s{1}{/;
				$flag=1;
			}
			elsif ($flag ==1) {
				if ($b==0) {
	        			if ($line =~ /lineage "/) {
						$line =~ /lineage "(.+)$/;
						$lineage = $1;
						if ($line =~ /.+"\s{1},/){ #whole entry on one line
							$line =~ /"(.+)"\s{1},/;
							$lineage = $1;
							print OUT "$lineage\t$gi\n";
							$flag=0;
							$a=0;
							$b=0;
						 	$c=1;
						}	
						else {
							print OUT "$lineage"; #first part of entry on one line
							$b = 1;
						}
					}
				}
				elsif ($b==1) {
					if ($line =~ /\s{1},/) { #last line of entry
						$line =~ /(.+)"\s{1},/;
						$lineage = $1;
						print OUT "$lineage\t$gi\n";
						$flag=0;
						$a = 0;
						$b = 0;
						$c=1;
					}
					else { #middle of entry
						print OUT "$line";
					}	
				}
			}
		}
		}
			
	}
	$c=0;
	close OUT;
	$i++;
}
