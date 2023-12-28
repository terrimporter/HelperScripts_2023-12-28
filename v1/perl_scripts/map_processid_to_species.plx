#!/usr/bin/perl
#May 30, 2012 by Terri Porter
#Script to associate processid from processid_NUC.fasta (from parse_Insecta.plx) with BIN_lineage.map (from parse_BOLD_html.plx)
#usage perl map_processid_to_species.plx processid_NUC.fasta BIN_lineage.map processid_BIN.map

use strict;
use warnings;

#declare fasta
my $i=0;
my $line;
my $j=0;
my $bin;
my $lineage;
my $line2;
my $species;
my $original;
my $new;
my $processid;
my $pid;
my $BIN;
my $scalar=0;
my $sp;
my $genus;
my $name;
my $binomial;

#declare array
my @bin;
my @pid;
my @lin;
my @line;
my @lineage;
my @species;
my @split;
my @binomial;

#declare hash
my %bin; #processid => bin
my %species; #bin => species1|species2 etc.

open (BIN, "<", $ARGV[2]) || die "Error cannot open processid_BIN.map: $!\n";
@bin = <BIN>;
close BIN;

#parse processid_BIN.map into hash
while ($bin[$i]) {
	$line = $bin[$i];
	chomp $line;
	
	@line = split(/\t/, $line);
	$processid = $line[0];
	$BIN = $line[1];

	$bin{$processid} = $BIN;
	$i++;

	$line=();
	@line=();
	$processid=();
	$BIN=();
}
$i=0;

open (LIN, "<", $ARGV[1]) || die "Error cannot open BIN_lineage.map: $!\n";
@lin = <LIN>;
close LIN;

#parse BIN_lineage.map to get just BIN and species rank associations (without the rest of the lineage)
while ($lin[$i]) {
	$line = $lin[$i];
	chomp $line;

	@line = split (/\t/, $line);
	$bin = $line[0];
	$lineage = $line[1];

	@lineage = split (/\|/, $lineage);
	while ($lineage[$j]) {
		$line2 = $lineage[$j];
		if ($line2 =~ /\[Species\]/) {
			$line2 =~ s/\[Species\]//;
			$species = $line2;

			@split = split(" ", $species);
			$genus = $split[0];

			#filter genus names to retain
			if ($genus =~ /^[a-z]/) {
				$genus=();
			}
			elsif ($genus =~ /\d+/) {
				$genus=();
			}
			
			$species = $split[1]; #ignore trinomial names
			#filter species names to retain

			if ($genus) {
				if ($species !~ /^[a-z]/) {
					print "$species\t";
					$species = ();
					print "\n";
					$name = $genus;
				}
					elsif ($species =~ /^\w{1,3}\d+/) {
					print "$species\t";
					$species =();
					$name = $genus;
					print "\n";
				}
				elsif ($species =~ /^\w{1,3}[A-Z]\S+/) {
					print "$species\t";
					$species = ();
					print "\n";
					$name = $genus;
				}
				elsif ($species =~ /sp\./) { #not identified to species rank so remove
					print "$species\t";
					if ($species =~ /sp\.\S+/) {
						$species =~ s/sp\.\S+//;
					}
					else {
						$species =~ s/sp\.//;
					}
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /_sp/) { #collapse all to just the species name without strain
					print "$species\t";
					if ($species =~ /_sp\d+/) {
						$species =~ s/_sp\d+//;
					}
					else {
						$species =~ s/_sp//;
					}
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /^sp[A-Z]/) {
					print "$species\t";
					$species = ();
					print "\n";
					$name = $genus;
				}
				elsif ($species =~ /^sp\d+/) {
					print "$species\t";
					$species =();
					print "\n";
					$name = $genus;
				}
				elsif ($species =~ /sp\d+\S*/) {
					print "$species\t";
					$species =~ s/sp\d+\S*//;
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /sp[A-Z]+\S*/) {
					print "$species\t";
					$species =~ s/sp[A-Z]+\S*//;
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /_cf/) { #not identified to species rank so remove
					print "$species\t";
					$species =~ s/_cf//;
					$species = ();
					$name = $genus;
					print "\n";
				}
				elsif ($species =~ /_\S+/) {
					print "$species\t";
					$species =~ s/_\S+//;
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /(biolep\S*|^\w{1,3}\d*$|species|geoBioLep\d+|^undet$)/) {
					print "$species\t";
					$species =();
					print "\n";
					$name = $genus;
				}
				elsif ($species =~ /[a-z]\d+/) {
					print "$species\t";
					if ($species =~ /[a-z]+[A-Z]+[a-z]+\d+/) {
						$species =~ /([a-z]+)[A-Z]+[a-z]+\d+/;
						$species = $1;
					}
					elsif ($species =~ /[a-z]\d+\S*/) {
						$species =~ s/\d+\S*//;
					}
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /[a-z][A-Z]/) {
					print "$species\t";
					if ($species =~ /[a-z][A-Z]+\S*/) {
						$species =~ s/[A-Z]+\S*//;
					}
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /\d+$/) { #exclude species char with numbered strains
					print "$species\t";
					if ($species =~ /[a-z]+[A-Z]+[a-z]*\d+$/) {
						$species =~ /([a-z]+)[A-Z]+[a-z]*\d+$/;
						$species = $1;
					}
					elsif ($species =~ /\d+/) {
						if ($species =~ /\-\d+/) {
							$species =~ s/\-\d+//g;
						}
						else {
							$species =~ s/\d//g;
						}
					}
					print "$species\n";
					$name = $genus." ".$species;
				}
				elsif ($species =~ /\./) {
					print "$species\t";
					if ($species =~ /^[a-z]{1,3}\.\S*/) {
						$species = ();
						$name = $genus;
					}
					elsif ($species =~ /^[a-z]+\-grp\.\S*/) {
						$species =~ /^([a-z]+)\-grp\.\S*/;
						$species = $1;
						$name = $genus." ".$species;
					}
					else {
						$species =~ /^([a-z]+)\.\S*/;
						$species = $1;
						$name = $genus." ".$species;
					}
					print "\n";
				}
				elsif ($species =~ /\|$/) {
					print "$species\t";
					$species =~ s/\|//g;
					print "$species\n";
					$name = $genus." ".$species;
				}
				else {
					$name = $genus." ".$species;
				}

		
				if (exists $species{$bin}) {
					$original = $species{$bin};
					$new = $original."|".$name;
					$species{$bin} = $new;
				}
				else {
					$species{$bin} = $name;
				}
			}
		}
		$j++;
		$line2=();
		$species=();
		$original=();
		$new=();
	}
	$j=0;

	$i++;
	$line=();
	@line=();
	$bin=();
	$lineage=();
	@lineage=();
	
}
$i=0;

open (PID, "<", $ARGV[0]) || die "Error cannot open processid_NUC.fasta: $!\n";
@pid = <PID>;
close PID;

open (OUT, ">>", "processid_species.map") || die "Error cannot open processid_species.map: $!\n";

#parse processid_NUC.fasta to get processid and species rank associations only (unix sort -u later to get a count of uniq species)
while ($pid[$i]) {
	$line = $pid[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$pid = $line;
		$BIN = $bin{$pid};
		if (exists $species{$BIN}) {
			$species = $species{$BIN};
			@species = split(/\|/, $species);
			$scalar = scalar(@species); #count number of species associated with BIN
		}

		if ($scalar == 1 ) {
			$binomial = $species[0];
			#print "binomial: $binomial\n";#test
			@binomial = split(" ",$binomial);
			$genus = $binomial[0];
			$sp = $binomial[1];
			if ($sp) { #only print if species epithet present
				print OUT "$pid\t$binomial\n";
			}
			$binomial=();
			@binomial=();
			$genus=();
			$sp=();
		}
		elsif ($scalar > 1) {

			while ($species[$j]) {
				$binomial = $species[$j];
				@binomial = split(" ",$binomial);
				$genus = $binomial[0];
				$sp = $binomial[1];
				if ($sp) { #only print if species epithet present
					print OUT "$pid\t$binomial\n";
				}
				$j++;
				$sp=();
				$binomial=();
				@binomial=();
				$genus=();
			}
			$j=0;

		}

	}

	$i++;
	$line=();
	$pid=();
	$BIN=();
	$species=();
	@species=();
	$scalar=0;
}
$i=0;
	
		


