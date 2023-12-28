#!/usr/bin/perl
#Dec. 22, 2009 by Terri Porter
#Script to turn newick trees from phylobayes or RaxML into a nexus formatted file that can be read by PAUP to create a consensus tree
#Usage $perl newicktonexus.plx < infile > outfile

use strict;
use warnings;

#declare variables
my $i=1; #counter
my $line;
my $val;
my $newval;

print "#NEXUS\nbegin trees;\n";
print "translate\n\t1 AA,\n\t2 AB,\n\t3 AM,\n\t4 AN,\n\t5 BB,\n\t6 BD,\n\t7 BE,\n\t8 CA,\n\t9 EP,\n\t10 GP,\n\t11 LB,\n\t12 MB,\n\t13 MC,\n\t14 MV,\n\t15 NE,\n\t16 PB,\n\t17 PE,\n\t18 RO,\n\t19 SC,\n\t20 SP,\n\t21 SR,\n\t22 UM;\n";

while (<STDIN>) {
	$line = $_;
	chomp($line);
	

	if (/;$/) {
		#first parse out two-letter taxon labels and replace with numbers
		$line =~ s/AA/1/;
		$line =~ s/AB/2/;
		$line =~ s/AM/3/;
		$line =~ s/AN/4/;
		$line =~ s/BB/5/;
		$line =~ s/BD/6/;
		$line =~ s/BE/7/;
		$line =~ s/CA/8/;
		$line =~ s/EP/9/;
		$line =~ s/GP/10/;
		$line =~ s/LB/11/;
		$line =~ s/MB/12/;
		$line =~ s/MC/13/;
		$line =~ s/MV/14/;
		$line =~ s/NE/15/;
		$line =~ s/PB/16/;
		$line =~ s/PE/17/;
		$line =~ s/RO/18/;
		$line =~ s/SC/19/;
		$line =~ s/SP/20/;
		$line =~ s/SR/21/;
		$line =~ s/UM/22/;

#		$line =~ /(\d+e-\d+)/; #change scientific notation to decimal
#		$val = $1;
#		$newval = sprintf("%.10f", $val);
#		$line;

#		$line =~ s/:\d+\.\d+e-\d+//g; ### also remove fields containing scientific notation ###
#		$line =~ s/:\d+\.\d+//g; ### remove branch lengths, since not used for bootstrapping ###
#		$line =~ s/:\d+e-\d+//g; ### whole numbers in scientific notation ###
#		$line =~ s/:\d+//g; ### check once more for whole numbers ###
		
		print "tree rep.",$i," = [&U] ",$line, "\n";
		$i++;
	}
	else {
		next;
	}
}

print "end;\n";
