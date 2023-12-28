#!/usr/bin/perl
#April 19, 2013 only need to remove '.1' from accessions so that TreeBase can auto sync with GenBank taxonomy...
#April 18, 2013 by Terri Porter
#Script to reformat the taxon labels for TreeBase submission
#usage perl fix_mesquite.plx file.nex

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $flag=0;
my $j=0;
my $label;
my $newlabel;
my $newline;
my $filename;
my $k;
my $previousline;
my $length;
my $space;
my $length2;
my $accession;

#declare array
my @in;
my @line;
my @label;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

#make changes
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^BEGIN\s{1}TAXA/) {
		$flag=1;
		$i++;
		next;
	}

	if ($flag==1) {
		if ($line =~ /TAXLABELS/) {
			$flag=2;
			$i++;
			next;
		}
	}

	if ($flag==2) {
		@line = split(' ', $line);

		while ($line[$j]) {
			$label = $line[$j];
			@label = split(/_/, $label);
			if ($label[0] !~ /\d+/) {
#				$newlabel = $label[0]."_".$label[1];
				$accession = $label[2];
				if ($accession =~ /^NC/) {
					$accession = $label[3];
					$accession =~ s/\.1//g;
					$newlabel = $label[0]."_".$label[1]."_".$label[2]."_".$accession;
				}
				else {
					$accession =~ s/\.1//g;
					$newlabel = $label[0]."_".$label[1]."_".$accession;
				}
			}
			elsif ($label[0] =~ /\d+/) {
				$newlabel = $label[0]."_".$label[1]."_".$label[2];
			}
			$line[$j] = $newlabel;
		
			$j++;
			$label=();
			@label=();
			$newlabel=();
		}
		$j=0;

		$newline = join(' ', @line);
		$newline = "\t\t".$newline;
		$in[$i] = $newline;
		
		$flag=3;
		$i++;
		next;
	}

	if ($flag == 3) {
		if ($line =~ /^BEGIN\s{1}CHARACTERS/) {
			$flag=4;
			$i++;
			next;
		}
	}
		
	if ($flag == 4) {
		if ($line =~ /MATRIX/) {
			$flag=5;
			$i++;
			next;
		}
	}

	if ($flag == 5 ) {
		if ($line !~ /;/) {
			if ($line ne "") {
				@line = split(/\s+/, $line);
				$label = $line[1]; #checked ok
#				print "label:$label\n";
				@label = split(/_/,$label);
				if ($label[0] !~ /\d+/) {
#					$newlabel = $label[0]."_".$label[1];
#					print "newlabel:$newlabel\n";
					$accession = $label[2];
#					print "accession:$accession\n";
					if ($accession =~ /^NC/) {
						$accession = $label[3];
						$accession =~ s/\.1//g;
						$newlabel = $label[0]."_".$label[1]."_".$label[2]."_".$accession;
					}
					else {
						$accession =~ s/\.1//g;
						$newlabel = $label[0]."_".$label[1]."_".$accession;
					}
				}
				elsif ($label[0] =~ /\d+/) {
					$newlabel = $label[0]."_".$label[1]."_".$label[2];
				}
					$line[1] = $newlabel;
					$length = length($newlabel);
					$length2 = 42-$length;
					$space = ' ' x $length2;
					$newline = '    '.$line[1].$space.$line[2];
#				print "newline:$newline\n";
				$newline = $newline;
				$in[$i] = $newline;
				$length=();
				$length2=();
				$space=();
				$i++;
				next;
			}
		}
		elsif ($line =~ /;/) {
			$flag=6;
			$i++;
			next;
		}
	}

	if ($flag == 6) {
		if ($line =~ /^BEGIN\s{1}TREES/) {
			$flag = 7;
			$i++;
			next;
		}
	}

	if ($flag == 7) {
		if ($line =~ /TRANSLATE/) {
			$flag = 8;
			$i++;
			next;
		}
	}

	if ($flag == 8) {
		if ($line !~ /TREE/) {
			@line = split(/\s+/,$line);
			$label = $line[2]; #check this
			@label = split(/_/,$label);
			if ($label[0] !~ /\d+/) {
#				$newlabel = $label[0]."_".$label[1].",";
				$accession = $label[2];
				if ($accession =~ /^NC/) {
					$accession = $label[3];
					$accession =~ s/\.1//g;
					$newlabel = $label[0]."_".$label[1]."_".$label[2]."_".$accession;
				}
				else {
					$accession =~ s/\.1//g;
					$newlabel = $label[0]."_".$label[1]."_".$accession;
				}
			}
			elsif ($label[0] =~ /\d+/) {
				$newlabel = $label[0]."_".$label[1]."_".$label[2].",";
			}
			$line[2] = $newlabel;
			$newline = join(' ',@line);
			$in[$i] = '       '.$newline;
			$i++;
			next;
		}
		else {
			$k = $i-1;
			$previousline = $in[$k];
			chomp $previousline;
			$previousline =~ s/,/;/g;
			$in[$k] = $previousline;
			$flag = 9;
			$i++;
			next;
		}
	}

	$i++;
}
$i=0;

$filename = $ARGV[0].".edited";

open (OUT, ">>", $filename) || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	print OUT "$line\n";

	$i++;
	$line=();
}
$i=0;

close OUT;
