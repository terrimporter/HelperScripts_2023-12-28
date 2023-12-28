#!/usr/bin/perl
#June 12, 2012 by Terri Porter
#Script to remove underscore from NC_number_number to NCnumber_number
#usage perl fix_NC.plx file.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $gb_id;

#declare array
my @in;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", "NC_fixed.fasta") || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		
		if ($line =~ /^NC_/) {
			$line =~ /^(NC)_(\w+)/;
			$gb_id = $1.$2;
		}
		else {
			$gb_id = $line;
		}
		print OUT ">$gb_id\n";
	}
	elsif ($line =~ /^\S/) {
		print OUT "$line\n";
	}
	
	$i++;
	$line=();
	$gb_id=();
}
$i=0;
close OUT;
