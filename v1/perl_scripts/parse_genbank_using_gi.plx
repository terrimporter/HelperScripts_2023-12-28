#!/usr/bin/perl
#January 27, 2010 by Terri Porter
#Updated Aug.31, 2010 to accept infile as an argument instead
#Script to parse genbank records retrieved using the script get_genbankrecord.plx and retrieve organism lineage.
#Usage #perl parse_genbankrecord.plx
#Script will automatically look for gb.outfile from get_genbankrecord.plx and will produce parsed_gb.outfile which is a tab-delimited file

use strict;
use warnings;

#declare variables
my $line;
my $a=0; #flags start of new record
my $gi;
my $b=0; #flags accession field
my $lineage;
my $c=0; #flags orgname field
my $d=0; #flags lineage field
my $k; #counter

#declare arrays
my @accessions;

open (IN, '<', $ARGV[0]) || die ("Error: $!\n");
@accessions = <IN>;
chomp(@accessions);

open (OUT, '>>', "parsed_genbankrecords.outfile") || die ("Error: $!\n");

foreach (@accessions) {
	$line = $_;
	chomp $line;
	if ($line =~ /Seq-entry ::=/) {
		$a = 1;
		$k++; #entry counter
		next;
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
			$b = 0;
			$c = 0;
		}
	}	
}
print $k; #test
