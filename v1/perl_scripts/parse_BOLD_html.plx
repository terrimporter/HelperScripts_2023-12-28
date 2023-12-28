#!/usr/bin/perl
#May 10, 2012 by Terri Porter
#Script to parse html files for each BIN (cluster ID) from BOLDSYSTEMS3
#usage perl parse_BOLD_html.plx

use strict;
use warnings;

#declare var
my $i=0;
my $filename;
my $BIN;
my $j=0;
my $line;
my $flag=0;
my $k;
my $line2;

#declare array
my @output;
my @filename;
my @in;

#declare hash
my %BIN;

@output = qx(ls | grep 'html\$');

#print "ArrayOutput: @output\n";#test

while ($output[$i]) { #remove any line endings from grepped filenames
	$filename = $output[$i];
	chomp $filename;
	#print "filename: $filename\n";#test

	@filename = split(/\./,$filename);
	$BIN = $filename[0];
	#print "BIN: $BIN\n";#test

	open (IN, "<", $filename) || die "Error cannot open html file: $!\n";
	@in = <IN>;
	#print "infile: @in\n";#test
	close IN;

	while ($in[$j]) { #parse individual html files
		$line = $in[$j];
		chomp $line;
		#print "$line\n";#test

		if ($line =~ /TAXONOMY/) {
			$flag=1;
			#print "flag:$flag\n";#test
		}
		elsif ($flag == 1) {
			if ($line =~ />Phylum</) {
				#print "found phylum\n";#test
				#print "$line\n";
				$k=$j+2;
				$line2 = $in[$k];
				chomp $line2;
				#print "$line2\n\n";
				if ($line2 =~ /\/table><table/) {#indicates multiple classifications
					process_multiple_classifications(\$line2,\$flag);
				}
				else { #process single classification
					process_single_classification(\$line2,\$flag);				
				}
				$flag=2;
			}
			$k=();
			$line2=();
		}
		elsif ($flag == 2) {
			if ($line =~ /> Class</) { #note space before Class
				#print "found Class\n";#test
				$k=$j+2;
				$line2 = $in[$k];
				chomp $line2;

				if ($line2 =~ /\/table><table/) {
					process_multiple_classifications(\$line2,\$flag);
				}
				else {
					process_single_classification(\$line2,\$flag);
				}
				$flag=3;
			}
			$k=();
			$line2=();
		}
		elsif ($flag == 3) {
			if ($line =~ />Order</) {
				$k=$j+2;
				$line2 = $in[$k];
				chomp $line2;
				if ($line2 =~ /\/table><table/) {
					process_multiple_classifications(\$line2,\$flag);
				}
				else {
					process_single_classification(\$line2,\$flag);
				}
				$flag=4;
			}
			$k=();
			$line2=();
		}
		elsif ($flag == 4) {
			if ($line =~ />Family</) {
				$k=$j+2;
				$line2 = $in[$k];
				chomp $line2;
				if ($line2 =~ /\/table><table/) {
					process_multiple_classifications(\$line2,\$flag);
				}
				else {
					process_single_classification(\$line2,\$flag);
				}
				$flag=5;
			}
			$k=();
			$line2=();
		}
		elsif ($flag == 5) {
			if ($line =~ />Genus</) {
				$k=$j+2;
				$line2 = $in[$k];
				chomp $line2;
				if ($line2 =~ /\/table><table/) {
					process_multiple_classifications(\$line2,\$flag);
				}
				else {
					process_single_classification(\$line2,\$flag);
				}
				$flag=6;
			}
			$k=();
			$line2=();
		}
		elsif ($flag == 6) {
			if ($line =~ />Species</) {
				$k=$j+2;
				$line2 = $in[$k];
				chomp $line2;
				if ($line2 =~ /\/table><table/) {
					process_multiple_classifications(\$line2,\$flag);
				}
				else {
					process_single_classification(\$line2,\$flag);
				}
				$flag=0;
			}
			$k=();
			$line2=();
		}
		$j++;
		$line=();
	}
	$j=0;
	$i++;
	$filename=();
	@filename=();
	$BIN=();
	@in=();
}
$i=0;

open (OUT, ">>", "BIN_lineage.map") || die "Error cannot open outfile: $!\n";

while (my ($BIN,$lineage) = each (%BIN)) {
	print OUT "$BIN\t$lineage\n";
}
close OUT;

####################

sub process_multiple_classifications {

#print "mult. classifiation needed\n";#test

#declare var
my $var = shift;#line2
my $var2 = $$var;
my $var3 = shift;#flag
my $var4 = $$var3;
my $tflag;
my $entry;
my $name;
my $oldvalue;
my $newvalue;

#declare array
my @var;

#translate flag from number to rank
if ($var4 == 1) {
	$tflag = "[Phylum]";
}
elsif ($var4 ==2) {
	$tflag = "[Class]";
}
elsif ($var4 == 3) {
	$tflag = "[Order]";
}
elsif ($var4 == 4) {
	$tflag = "[Family]";
}
elsif ($var4 == 5) {
	$tflag = "[Genus]";
}
elsif ($var4 == 6) {
	$tflag = "[Species]";
}

@var = split(/NOWRAP>/, $var2);

foreach $entry (@var) {
	if ($entry =~ /^\s+<table/) {
		next;
	}
	else {
		#print "entry:$entry\n";
		$entry =~ /(.+)&nbsp;/;
		$name = $1;

		if (exists $BIN{$BIN}) {
			$oldvalue = $BIN{$BIN};
			$newvalue = $oldvalue."|".$name.$tflag;
			$BIN{$BIN} = $newvalue;
		}
		else {
			$BIN{$BIN} = $name.$tflag;
		}
	}
}

}

####################

sub process_single_classification {

#declare var
my $var = shift; #grab first element of array passed as reference into sub
my $var2 = $$var;
my $var3 = shift; #flag
my $var4 = $$var3;
my $tflag;
my $name;
my $oldvalue;
my $newvalue;

#translate flag from number to rank
if ($var4 == 1) {
	$tflag = "[Phylum]";
}
elsif ($var4 ==2) {
	$tflag = "[Class]";
}
elsif ($var4 == 3) {
	$tflag = "[Order]";
}
elsif ($var4 == 4) {
	$tflag = "[Family]";
}
elsif ($var4 == 5) {
	$tflag = "[Genus]";
}
elsif ($var4 == 6) {
	$tflag = "[Species]";
}

#print "var2:$var2\n";#test
if ($var2 =~ /NOWRAP>(\w+|\w+\s+\w+|\w+\s+\w+\s+\w+)&nbsp;/) {
	$var2 =~ /NOWRAP>(\w+|\w+\s+\w+|\w+\s+\w+\s+\w+)&nbsp;/;
	$name = $1;
	#print "name:$name\n";#test

	if (exists $BIN{$BIN}) {
		$oldvalue = $BIN{$BIN};
		$newvalue = $oldvalue."|".$name.$tflag;
		$BIN{$BIN} = $newvalue;
	}	
	else {
		$BIN{$BIN} = $name.$tflag;
	}
}

}
