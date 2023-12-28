#!/usr/bin/perl
#Oct.5, 2010 by Terri Porter
#Script to produce tabular output from clone.html files instead of the index.html file, concentrating on species level hits
#usage $ perl parse_clone_html_files.plx

use strict;
use warnings;

#declare var
my $file;
my $line;
my $id;
my $i=0;
my $j=0;
my $flag=0;
my $genus;
my $species;
my $pp;

#declare array
my @files;
my @current_file;
my @line;

my $dir = "/home/terri/SAP_data/project-14.47-10.04.10_Kurtzman_fillinall_minidentity/html/clones/";
opendir (DH, $dir) || die "Error:$!\n";
@files = readdir (DH);
#print "@files\n";#test

#open (OUT,"<","clone_summary.txt") || die ("Error:$!\n");

while ($files[$i]) {
	$file = $files[$i];
	if ($file =~ /Saccharomyces/) {

		open (FH, $file) || die ("Error:$!\n");
		@current_file = <FH>;
		while ($current_file[$j]){
			$line = $current_file[$j];
			chomp $line;
			if ($flag==0 && $line =~ /<h1>Sequence:/) {
				$line =~ /<h1>Sequence: (\S+)<\/h1>/;
				$id = $1;
				print "$id\t";
				$flag=1;
			}
			elsif ($flag==1) {
			       if ($line =~ /species/) {
				$flag=2;
				$line =~/\+(\S+)\s{1}(\S+)\s{1}\S{1}species\S{1}\s{1}(\S+)%/;
				$genus = $1;
				$species = $2;
				$pp = $3;
				print "$genus\t$species\t$pp\n";
				}
			}
			elsif ($flag==2 && $line =~ /species/) {
				#print "$id\t";
				$line =~ /\+(\S+)\s{1}(\S+)\s{1}\S{1}species\S{1}\s{1}(\S+)%/;
				$genus = $1;
				$species = $2;
				$pp = $3;
				print "$id\t$genus\t$species\t$pp\n";
			}
			elsif ($flag==2 && $line =~ /<\/pre>/) {
				$flag=0;
			}
			$j++;
		}
		$j=0;
		$i++;
		close FH;
	}
	else {
		$i++;
		next;
	}
}


