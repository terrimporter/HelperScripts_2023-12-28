#!/usr/bin/perl
#Jan. 17, 2012 edited to parse NBC output from locally installed fungal LSU classifier and filter results by statistical cutoff (0.80 for 250 bp +; otherwise 0.50)
#Sept.30, 2011 edited to parse NBC output and filter results by statistical cutoff (80 default)
#new usage perl get_recovery2.plx features.txt.parsed NBC_download.fasta
#Sept.29, 2011 by Terri Porter
#Script to compare features.txt.parsed from test set with megan_hit_lineage.txt.mapped
#usage perl get_recovery.plx features.txt.parsed gb.query gi.query megan_hit_lineage.txt.mapped

use strict;
use warnings;

#declare var
my $CUTOFF = 0.8;  #####SET CUTOFF HERE##### use 0.8 for sequences 250 bp +, otherwise 50 for short seqs
#declare array
my @ref;
my @gb;
my @gi;
my @nbc;

#declare hash
my %map_gb_gi;#this hash keyed with gb
my %ref_organism;#each hash below keyed with gb
my %ref_kingdom;
my %ref_phylum;
my %ref_class;
my %ref_order;
my %ref_family;
my %ref_genus;
my %ref_species;
my %kingdom;
my %phylum;
my %class;
my %order;
my %family;
my %genus;
my %species;
my %kingdom_stat;
my %phylum_stat;
my %class_stat;
my %order_stat;
my %family_stat;
my %genus_stat;
my %species_stat;

open (REF,"<",$ARGV[0]) || die "Error cannot open features.txt.parsed:$!\n";
@ref = <REF>;
close REF;

open (NBC,"<",$ARGV[1]) || die "Error cannot open nbc_output:$!\n";
@nbc = <NBC>;
close NBC;

parse_ref_lineage();

parse_NBC_lineage();

compare_classifications();

##########parse lineage info from test set
sub parse_ref_lineage {

#declare var
my $b=0;
my $line;
my $ref_gb;
my $ref_gi;
my $ref_organism;
my $ref_strain;
my $ref_isolate;
my $ref_length;
my $ref_classification;
my $ref_genus;
my $ref_species;
my $x;
my $ref_kingdom;
my $ref_phylum;
my $ref_class;
my $ref_order;
my $ref_family;

#declare array
my @line;
my @organism;
my @classification;

while ($ref[$b]) {
	$line = $ref[$b];
	chomp $line;
	@line =	split(/\t/,$line);
	$ref_gb = $line[0];
	#$ref_gi = $map_gb_gi{$ref_gb};
	
	$ref_organism = $line[1];
	$ref_organism{$ref_gb}= $ref_organism;
	@organism = split(/ /,$ref_organism);
	$ref_genus = $organism[0];
	$ref_genus{$ref_gb} = $ref_genus;
	$ref_species = $organism[1];
	$ref_species{$ref_gb} = $ref_species;
	
	$ref_strain = $line[2];
	$ref_isolate = $line[3];
	$ref_length = $line[4];
	$ref_classification = $line[5];
	@classification = split(/ /,$ref_classification);

	foreach $x (@classification) {
		if ($x eq "Fungi") {
			$ref_kingdom = $x;
			$ref_kingdom{$ref_gb} = $ref_kingdom;
		}
		elsif ($x =~ /mycota$/) {
			$ref_phylum = $x;
			$ref_phylum{$ref_gb} = $ref_phylum;
		}
		elsif ($x =~ /etes$/) {
			$ref_class = $x;
			$ref_class{$ref_gb} = $ref_class;
		}
		elsif ($x =~ /ales$/) {
			$ref_order = $x;
			$ref_order{$ref_gb} = $ref_order;
		}
		elsif ($x =~ /aceae$/) {
			$ref_family = $x;
			$ref_family{$ref_gb} = $ref_family;
		}
	}
	$b++;
	@line=();
	@organism=();
	@classification=();
	$ref_gb=();
	$ref_gb=();
	$ref_organism=();
	$ref_genus=();
	$ref_species=();
	$ref_strain=();
	$ref_isolate=();
	$ref_length=();
	$ref_classification=();
	$ref_kingdom=();
	$ref_phylum=();
	$ref_class=();
	$ref_order=();
	$ref_family=();
	
}

}

##########parse lineage info from seqs classified by megan via taxonids and [r] taxid2names

sub parse_NBC_lineage {

#declare var
my $c=0;
my $line;
my $flag=0;
my $fragment1;
my $ref_gb;
my $d=0;
my $x;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $e;
my $species;
my $genus;
my $stat;
my $prob;
my $f="nil";
my $g="nil";

#declare array
my @line;
my @fragment1;

while ($nbc[$c]) {
	$line = $nbc[$c];
	chomp $line;
	#print "nbc_line:$line\n";#test
	if ($line =~ /^\w+/) { #not ^Details:
		$flag=1;
	}
	if ($flag ==1) {
		@line = split(/\t/,$line); #not ; delimited
		#$ref_gb = shift @line;
		$ref_gb = $line[0];
		print "gb: $ref_gb\t";
	
		#while ($line[$d]) {
		#	$x = $line[$d];
		#	$x =~ s/^\s{1}//;
		#	if ($x =~ /&amp;apos;/) {
		#		$x =~ s/&amp;apos;//;
		#	}
			$kingdom = $line[5]; #domain id
		#	if ($x eq "Fungi") {
		#		$kingdom = $x;
				$kingdom{$ref_gb} = $kingdom;
		#		$e = $d+1;
		#		$stat = $line[7];
		#		$stat =~ /^\s{1}(\d+)%/;
				$prob = $line[7];
				$kingdom_stat{$ref_gb} = $prob;
#				print "kingdom: $kingdom\t$prob\t";
		#		$e=();
		#		$stat=();
				$prob=();
		#	}
		#	elsif ($x =~ /mycota$/) {
		#		if ($x !~ /(mitosporic|Euglenida)/) {
					$phylum =$line[8];
					$phylum{$ref_gb} = $phylum;
		#			$e = $d+1;
		#			$stat = $line[$e];
		#			$stat =~ /^\s{1}(\d+)%/;
					$prob = $line[10];
					$phylum_stat{$ref_gb} = $prob;
#					print "phylum: $phylum\t$prob\t";
		#			$e=();
		#			$stat=();
					$prob=();
		#		}
		#	}
		#	elsif ($x =~ /etes$/) {
				$class =$line[11];
				$class{$ref_gb} = $class;
#				$e = $d+1;
#				$stat = $line[$e];
#				$stat =~ /^\s{1}(\d+)%/;
				$prob = $line[13];
				$class_stat{$ref_gb} = $prob;
#				print "class: $class\t$prob\t";
#				$e=();
#				$stat=();
				$prob=();
#			}
#			elsif ($x =~ /ales$/) {
				$order =$line[14];
#				if ($x !~ /(mitosporic|Choanoflagellida)/) {
					$order{$ref_gb} = $order;
#					$e = $d+1;
#					$stat = $line[$e];
#					$stat =~ /^\s{1}(\d+)%/;
					$prob = $line[16];
					$order_stat{$ref_gb} = $prob;
#					print "order: $order\t$prob\t";
#					$e=();
#					$stat=();
					$prob=();
#				}
#			}
#			elsif ($x =~ /aceae$/) {
				$family =$line[17];
#				if ($x !~ /(mitosporic|Codonosigidae)/) {
					$family{$ref_gb} = $family;
#					$e = $d+1;
#					$stat = $line[$e];
#					$stat =~ /^\s{1}(\d+)%/;
					$prob = $line[19];
					$family_stat{$ref_gb} = $prob;
#					print "famiy: $family\t$prob\t";
#					$stat=();
					$prob=();
#					$f=$e+1;#genus
#					$g=$e+2;#genus_stat
#					$e=();
#				}
#			}
#			elsif ($f ne "nil") {
				$genus = $line[20];
#				$genus =~ s/^\s{1}//;
#				if ($genus !~ /(incertae sedis|complex|group)/) {
					$genus{$ref_gb} = $genus;
#					if ($g ne "nil") {
#						$stat = $line[$g];
#						$stat =~ /^\s+(\d+)%/;
						$prob = $line[22];
						$genus_stat{$ref_gb} = $prob;
#						print "genus: $genus\t$prob\n";
#						$e=();
#						$f="nil";
#						$g="nil";
#						$stat=();
						$prob=();
#					}
#				}
#			}
#			$d++;
#		}
#		$d=0;
#		$e=();
#		$f="nil";
#		$g="nil";
		@line=();
		$kingdom=();
		$phylum=();
		$class=();
		$order=();
		$family=();
		$species=();
		$genus=();
	}
	$c++;
}	

}

##########For each gi, and for each rank, compare classifiation hashes to calculate recovery and erroneous recovery
sub compare_classifications {

#declare var
my $gb;
my $ref_kingdom;
my $kingdom;
my $coverage=0;
my $correct=0;
my $incorrect=0;
my $ref_phylum;
my $phylum;
my $ref_class;
my $class;
my $ref_order;
my $order;
my $ref_family;
my $family;
my $ref_genus;
my $genus;
my $ref_species;
my $species;
my $kingdom_stat;
my $phylum_stat;
my $class_stat;
my $order_stat;
my $family_stat;
my $genus_stat;

#declare array

open (OUT,">>","comparison.txt") || die "Error cannot write to comparison.txt:$!\n";

while (($gb, $kingdom) = each(%kingdom) ) {	
	$kingdom = $kingdom{$gb};
	$ref_kingdom = $ref_kingdom{$gb};
	if ($kingdom =~ /^\S+/) {
		$coverage++;
		if ($ref_kingdom eq $kingdom) {
			$kingdom_stat = $kingdom_stat{$gb};
			print "ref kingdom:$ref_kingdom\tkingdom:$kingdom\tkingdom_stat:$kingdom_stat\n";
			if ($kingdom_stat >= $CUTOFF) {
				$correct++;
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_kingdom ne $kingdom) {
			if ($kingdom_stat >= $CUTOFF) {
				$incorrect++;
			}
			else {
				$coverage--;
			}
		}
	}
}

print OUT "rank\tcorrect\tincorrect\tcoverage\n";
print OUT "kingdom\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

while (($gb, $phylum) = each(%phylum) ) {	
	$ref_phylum = $ref_phylum{$gb};
	$phylum = $phylum{$gb};
	if ($phylum =~ /^\S+/) {
		$coverage++;
		if ($ref_phylum eq $phylum) {
			$phylum_stat = $phylum_stat{$gb};
			if ($phylum_stat >= $CUTOFF) {
				$correct++;
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_phylum ne $phylum) {
			if ($phylum_stat >= $CUTOFF ) {
				$incorrect++;
			}
			else {
				$coverage--;
			}
		}
	}
}

print OUT "phylum\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

while (($gb, $class) = each(%class) ) {	
	$ref_class = $ref_class{$gb};
	$class = $class{$gb};
	if ($class =~ /^\S+/) {
		$coverage++;
		if ($ref_class eq $class) {
			$class_stat = $class_stat{$gb};
			if ($class_stat >= $CUTOFF) {
				$correct++;
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_class ne $class) {
			if ($class_stat >= $CUTOFF ) {
				$incorrect++;
			}
			else {
				$coverage--;
			}
		}
	}
}

print OUT "class\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

while (($gb, $order) = each(%order) ) {	
	$ref_order = $ref_order{$gb};
	$order = $order{$gb};
	if ($order =~ /^\S+/) {
		$coverage++;
		if ($ref_order eq $order) {
			$order_stat = $order_stat{$gb};
			if ($order_stat >= $CUTOFF) {
				$correct++;
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_order ne $order) {
			if ($order_stat >= $CUTOFF) {
				$incorrect++;
			}
			else {
				$coverage--;
			}
		}
	}
}

print OUT "order\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

while (($gb, $family) = each(%family) ) {	
	$ref_family = $ref_family{$gb};
	$family = $family{$gb};
	if ($family =~ /^\S+/) {
		$coverage++;
		if ($ref_family eq $family) {
			$family_stat = $family_stat{$gb};
			if ($family_stat >= $CUTOFF) {
				$correct++;
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_family ne $family) {
			if ($family_stat >= $CUTOFF ) {
				$incorrect++;
			}
			else {
				$coverage--;
			}
		}
	}
}

print OUT "family\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

while (($gb, $genus) = each(%genus) ) {	
	$ref_genus = $ref_genus{$gb};
	$genus = $genus{$gb};
	#print "Ref genus: $ref_genus\tgenus: $genus\n";#test
	if ($genus =~ /^\S+/) {
		$coverage++;
		if ($ref_genus eq $genus) {
			$genus_stat = $genus_stat{$gb};
			if ($genus_stat >= $CUTOFF) {
				$correct++;
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_genus ne $genus) {
			if ($genus_stat >= $CUTOFF) {
				$incorrect++;
			}
			else {
				$coverage--;
			}
		}
	}
}

print OUT "genus\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

close OUT;

}
