#!/usr/bin/perl
#January 27, 2010 by Terri Porter
#Updated Sept.13, 2010 to also retrieve species from binomial field
#Updated Aug.31, 2010 to accept infile as an argument instead
#Script to parse genbank records retrieved using the script get_genbankrecord.plx and retrieve organism lineage.
#Usage #perl parse_genbankrecord.plx
#Script will automatically look for gb.outfile from get_genbankrecord.plx and will produce parsed_gb.outfile which is a tab-delimited file

use strict;
use warnings;

#declare variables
my $line;
my $a=0; #flags start of new record
my $e=0;
my $species;
my $gi;
my $b=0; #flags accession field
my $lineage;
my $c=0; #flags orgname field
my $d=0; #flags lineage field
my $k; #counter
my $z; #index

#declare arrays
my @accessions;
my @line;
my @lineage;

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
@accessions = <IN>;
chomp(@accessions);

open (OUT, '>>', "temp.outfile") || die ("Error: $!\n");

foreach (@accessions) {
	$line = $_;
	chomp $line;
	if ($line =~ /Seq-entry ::=/) {
		$a = 1;
		$k++; #entry counter
		next;
	}
	elsif ($a==1 && $e==0) {
		if ($line =~ /species "/){
			$line =~ /species "(\w+)"/;
			$species = $1;
			print OUT "$species\t";
			$e=1;
		}
	}
	elsif ($a == 1 && $b==0) {
	        if ($line =~ /lineage "/) {
			$line =~ /lineage "(.+)$/;
			$lineage = $1;
			if ($line =~ /.+"\s{1},/){ #whole entry on one line
				$line =~ /"(.+)"\s{1},/;
				$lineage = $1;
				print OUT "$lineage\t";
				$b = 0;
				$c = 1;
			}	
			else {
				print OUT "$lineage"; #first part of entry on one line
				$b = 1;
			}
		}
	}
	elsif ($b == 1 && $c==0) {
		if ($line =~ /\s{1},/) { #last line of entry
			$line =~ /(.+)"\s{1},/;
			$lineage = $1;
			print OUT "$lineage\t";
			$c = 1;
		}
		else { #middle of entry
			print OUT "$line";
		}	
	}
	elsif ($c == 1 && $a==1){
	       if ($line =~ /gi\s{1}\d+\s{1}\}\s{1},/) {
			$line =~ /gi\s{1}(\d+)\s{1}\}\s{1},/;
			$gi = $1;
			print OUT "$gi\n";
			$a = 0;
			$e = 0;
			$b = 0;
			$c = 0;
		}
	}	
}
print $k; #test
close IN;
close OUT;

open (IN2,"<","temp.outfile") || die ("Error:$!\n");
open (OUT2,">>","parsed_genbankrecords_with_species.outfile") || die ("Error:$!\n");

while (<IN2>){
	$line = $_;
	chomp $line;
	@line = split(/\t/,$line);
	$gi = $line[2];
	print OUT2 "$gi\t";
	$lineage = $line[1];
	@lineage = split (/;/,$lineage);
	foreach $z (@lineage) {
		print OUT2 "\t$z\t";
	}
	$species = $line[0];
	print OUT2 "$species\n";
}
close IN2;
close OUT2;
unlink("temp.outfile");

