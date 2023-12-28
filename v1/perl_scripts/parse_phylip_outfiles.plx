#!/usr/bin/perl
#Edited Nov. 6, 2012 to work right
#Dec. 8, 2011 by Terri Porter
#Script to parse outfiles from run_phylip.sh
#usage perl parse_phylip_outfiles.plx

use strict;
use warnings;

#declare var
my $dir;
my $i=0;
my $file;
my $id;
my $path_to_infile;
my $j=0;
my $line;
my $flag=0;
my $num;
my $taxon;
my $tot_bootstrap;
my $pattern1;
my $pattern2;
my $pattern;
my $k=0;
my $symbol;
my $l;
my $support;
my $m=0;
my $ref_taxon;
my $path_to_outfile;

#declare array
my @files;
my @infile;
my @taxa;
my @pattern;
my @reference_taxa = ("LB","UM","SR","AN","SC","AA","MV","EP","RO","MC","PB","NE","SP","AB","BD","PE","GP","CA","BB","BE","AM","MB");

#declare hash
my %taxa;
my %pattern;
my %taxa2;
my %pattern2;

print "Enter path to directory containing orthoid.outfiles including final / : $!\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error cannot open directory: $!\n";
@files = readdir (DIR);
closedir (DIR);

$path_to_outfile = $dir."ortho.clades";
open (OUT, ">>", $path_to_outfile) || die "Error cannot write to outfile: $!\n";
print OUT "OrthoID\tClade\tBootstrap support\tTotal bootstraps\n";

while ($files[$i]) {
	$file = $files[$i];
	#print "checking outfiles\n";
	if ($file =~ /\d+\.outfile/)  {
		$file =~ /(\d+)\.outfile/;
		#print "$file\n";
		$id = $1;
		#print "$id\t";
		$path_to_infile = $dir.$file;

		open (IN, "<", $path_to_infile) || die "Error cannot open outfile: $!\n";
		@infile = <IN>;
		close IN;

		while ($infile[$j]) {
			$line = $infile[$j];
			#print "parsing infile\n";
			if ($flag==0 && $line =~ /^Species in order:/) {
				$flag=1;
				#print "found species line\n";
			}
			
			elsif ($flag==1 && $line =~ /\d+\.\s+\w+/) {
				$line =~ /(\d+)\.\s+(\w+)/;
				$num = $1;
				$taxon = $2;
				#print "num: $num\t taxon: $taxon\n";
				push(@taxa, $taxon); #just to easily get total number of taxa???
				$taxa{$num} = $taxon; #otherise use the hash
				$taxa2{$taxon} = $num; #to test for presence
			}
			elsif ($flag==1 && $line =~ /^\s+$/g) {
				$num = "blank";
				$taxon = "blank";
				print "num: $num\t taxon: $taxon\n";
			}
			elsif ($flag==1 && $line =~ /Sets included in the consensus tree/) {
				$flag=2;
				#print "found sets line\n"
			}
			elsif ($flag==2 && $line =~ /How many times out of\s+\d+/) {#fixed regex
				$line =~ /How many times out of\s+(\d+)/;#fixed regex
				$tot_bootstrap = $1;
				print "tot bootstrap: $tot_bootstrap\n";
			}
			elsif ($flag==2 && $line =~ /^(\*|\.){1,10}\s{0,1}/) {
				$line =~ /^((\*|\.){1,10})\s{0,1}/;
				$pattern1 = $1;
				$line =~ s/^(\*|\.){1,10}\s{0,1}//;
				if ($line =~ /(\*|\.){0,10}\s+/) {
					$line =~ /((\*|\.){0,10})\s+/;
					$pattern2 = $1;
					$line =~ s/^(\*|\.){0,10}\s+//;#added ^
					if ($line =~ /\d+\.\d+/) {
						$line =~ /(\d+)\.\d+/;
						$support = $1;
					}
				}
				print "pattern1: $pattern1\tpattern2: $pattern2\tsupport: $support\ttot_bootstrap: $tot_bootstrap\n";
				$pattern=$pattern1.$pattern2;
				@pattern = split(//,$pattern);
				print "pattern array: @pattern\t support: $support\n";
				while ($pattern[$k]) {
						$symbol = $pattern[$k];
						$l = $k+1; #to match taxon number from above
						$pattern{$l} = $symbol;
						$pattern2{$symbol} = $l;
						print "num: $l\tsymbol: $symbol\n";
						$k++;
						$l=(); #added
				}
				$k=0;
				print OUT "$id\t";
				#reformat patterns to be in same species order (*, ., -)
				while ($reference_taxa[$m]) {
					$ref_taxon = $reference_taxa[$m];
					if (exists $taxa2{$ref_taxon}) { #added exists
						$num = $taxa2{$ref_taxon};
#						if (exists $pattern{$num}) { #added exists
							$symbol = $pattern{$num};
							print OUT $symbol;
#						}
#						else { #added
#							print OUT "-";
#						}
					}
					else {
						print OUT "-";
					}
					$m++;
				}
				$m=0;
				print OUT "\t$support\t$tot_bootstrap\n";
			}
			elsif ($flag==2 && $line =~ /Sets NOT included in consensus tree/) {
				$flag=3;
			}
			$j++;
			@pattern=();
			$pattern1=();
			$pattern2=();
			$support=();
		}
		$flag=0;
		$j=0;
		$num=();
		$taxon=();
		@taxa=();
		%taxa=();
		%taxa2=();
		$tot_bootstrap=();
		$pattern=();
		%pattern=();
		%pattern2=();
	}
	$i++;
}
$i=0;
