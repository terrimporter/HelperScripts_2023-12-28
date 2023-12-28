#!/usr/bin/perl
#Sept. 7, 2011 by Terri Porter
#Script to parse web of science abstract for ITS and LSU keywords
#usage perl parse_abstract.plx ID_ABS.txt

#declare var
my $i=0;
my $line;
my $id;
my $abstract;

#declare array
my @in;

use strict;
use warnings;

open (IN,"<",$ARGV[0]) || die ("Error cannot open ID_ABS.txt:$!\n");
@in = <IN>;

open (OUT,">>","ITS_LSU.txt") || die ("Error cannot write to outfile: $!\n");

print OUT "ID\tITS\tLSU\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^\d+/) {
		if ($line =~ /^\d+/) {
			$line =~ /^(\d+)/;
			$id = $1;
			print OUT "$id\t";
			$line =~ /^\d+\t(.+)/;
			$abstract = $1;
			#print "$abstract\n";
			if ($abstract =~ /(ITS|internal transcribed)/) {
				print OUT "1\t";
				if ($abstract =~ /(LSU|lsu |lsu-|28S|28s|26S|26s|25S|25s|large subunit|large ribosomal)/) {
					print OUT "1\n";
				}
				else {
					print OUT "0\n";
				}
			}
			elsif ($abstract =~ /(LSU|lsu |lsu-|28S|28s|26S|26s|25S|25s|large subunit|large ribosomal)/) {
				print OUT "0\t1\n";
			}
			else {
				print OUT "0\t0\n";
			}
		}
		elsif ($line =~ /^\d+$/) {
			$line =~ /^(\d+)$/;
			$id = $1;
			print OUT "$id\n";	
		}
	}

	$i++;
	$abstract="";
}
