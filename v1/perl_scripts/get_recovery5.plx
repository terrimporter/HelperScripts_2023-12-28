#!/usr/bin/perl
#March 5, 2012 edit to read megan_hit_lineage.txt.mapped file that contains gb instead of gi properly ##originally downloaded NBC classifications, LCA parsed them in MEGAN, grabbed taxonids and lineages using [r], now get stats (correct, incorrect, coverage)
#Sept.29, 2011 by Terri Porter
#Script to compare features.txt.parsed from test set with megan_hit_lineage.txt.mapped
#usage perl get_recovery.plx features.txt.parsed gb.query gi.query megan_hit_lineage.txt.mapped

use strict;
use warnings;

#declare var

#declare array
my @ref;
my @gb;
my @gi;
my @megan;

#declare hash
my %map_gb_gi;#this hash keyed with gb
my %ref_organism;#each hash below keyed with gi
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

open (REF,"<",$ARGV[0]) || die "Error cannot open features.txt.parsed:$!\n";
@ref = <REF>;
close REF;

open (GB,"<",$ARGV[1]) || die "Error cannot open gb.query:$!\n";
@gb=<GB>;
close GB;

open (GI,"<",$ARGV[2]) || die "Error cannot open gi.query:$!\n";
@gi=<GI>;
close GI;

open (MEGAN,"<",$ARGV[3]) || die "Error cannot open megan_hit_lineage.txt.mapped:$!\n";
@megan = <MEGAN>;
close MEGAN;

make_map();

parse_ref_lineage();

parse_megan_lineage();

compare_classifications();

##########create hash to map gb and gi from test set
sub make_map {

#declare var
my $a=0;
my $gb;
my $gi;

while($gb[$a]) {
	$gb = $gb[$a];
	chomp $gb;
	$gb =~ /(\w+)\.\d+/;
	$gb = $1;
	$gi = $gi[$a];
	chomp $gi;
	$map_gb_gi{$gb} = $gi;
	$a++;
}

}

#while(my($key,$value) = each(%map_gb_gi)) {
#	print "key:$key\tvalue:$value\n";
#}#test

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
#	$ref_gi = $map_gb_gi{$ref_gb}; ##remove this line because want to index by $ref_gb
	
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
	$ref_gi=();
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

#while(my($key,$value) = each(%ref_species)) {
#	print "key:$key\tvalue:$value\n";
#}#test

}

##########parse lineage info from seqs classified by megan via taxonids and [r] taxid2names

sub parse_megan_lineage {

#declare var
my $c=0;
my $line;
my $fragment1;
my $ref_gi;
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

#declare array
my @line;
my @fragment1;

while ($megan[$c]) {
	$line = $megan[$c];
	chomp $line;
	#$line =~ s/"$//;
	#print "megan_line:$line\n";#test
	
	@line = split(/;/,$line);
	$fragment1 = shift @line;
	@fragment1 = split(/\s+/,$fragment1);
	$ref_gb = $fragment1[0];##changed from $ref_gi to $ref_gb
	
	while ($line[$d]) {
		$x = $line[$d];
		$x =~ s/^\s{1}//;
		$x =~ s/^&amp;//g;
		$x =~ s/^apos;//g;
		
		if ($x =~ "Fungi") {
			if ($x =~ /Fungi$/) {
				$kingdom = $x;
				$kingdom{$ref_gb} = $kingdom;
			}
			elsif ($x =~ /Fungi\./) {
				$x =~ /(Fungi)\.$/;
				$kingdom = $1;
				$kingdom{$ref_gb} = $kingdom;
				$d++;
				next;
			}
			#print "kingdom: $kingdom\n";#test
		}
		elsif ($x =~ /mycota/) {
			if ($x !~ /(mitosporic|Euglenida)/) {
				if ($x =~ /mycota$/) {
					$phylum =$x;
					$phylum{$ref_gb} = $phylum;
				}
				elsif ($x =~ /mycota\.$/) {
					$x =~ s/\.$//;
					$phylum = $x;
					$phylum{$ref_gb} =$phylum;
					$d++;
					next;
				}
				#print "phylum: $phylum\n";#test
			}
		}
		elsif ($x =~ /etes/) {
			if ($x !~ /mitosporic/) {
				if ($x =~ /etes$/) {
					$class = $x;
					$class{$ref_gb} = $class;
				}
				elsif ($x =~ /etes\.$/) {
					$x =~ s/\.$//;
					$class =$x;
					$class{$ref_gb} = $class;
					$d++;
					next;
				}
				#print "class: $class\n";#test
			}
		}
		elsif ($x =~ /ales/) {
			if ($x !~ /(mitosporic|Choanoflagellida)/) {
				if ($x =~ /ales$/) {
					$order = $x;
					$order{$ref_gb} = $order;
				}
				elsif ($x =~ /(ales)\.$/) {
					$x =~ s/\.$//;
					$order =$x;
					$order{$ref_gb} = $order;
					$d++;
					next;
				}
				#print "order: $order\n";
			}
		}
		elsif ($x =~ /aceae/) {
			if ($x !~ /(mitosporic|Codonosigidae)/) {
				if ($x =~ /aceae$/) {
					$family = $x;
					$family{$ref_gb} = $family;
				}
				elsif ($x =~ /aceae\.$/) {
					$x =~ s/\.$//;
					$family =$x;
					$family{$ref_gb} = $family;
					$d++;
					$genus =();
					$genus{$ref_gb} = $genus;
#					$species = ();
#					$species{$ref_gb} = $species;
					next;
				}
				#print "family: $family\n";
				#$e=$d+1;
				#$genus = $line[$e];
				#$genus =~ s/^\s{1}//;
				#if ($genus !~ /(incertae sedis|complex|group|mitosporic)/) {
					#print "genus: $genus\n";#test
					#if ($genus =~ /\w+$/) {
					#	$genus{$ref_gi} = $genus;
					#}
				#	if ($genus =~ /\w+\.$/) {
				#		$genus =~ s/\.$//;
				#		$genus{$ref_gi} = $genus;
				#		print "genus:$genus\n";#test
				#		$d++;
				#		next;
				#	}
				#}
			}
		}
#		elsif ($x =~ /\.\.$/) {
#			$x =~ s/^\.+//;
#			$x =~s/^\s+//;
#			$x =~ s/\.+$//;
#			if ($x =! /(Fungi|mycota$|etes$|ales$|aceae$|mitosporic|^[a-z])/) { #make sure that starts with capital letter if it's a genus; also make sure that it's not another rank that I can check for....
#				$species = $x;
#				$species{$ref_gb} = $species;
#			}
			#print "species:$species\n";
#		}
		else {
			if ($x !~ /\.$/) {
				$x =~ s/^\s+//;
#				$x =~ s/\.+//;
				if ($x !~ /(incertae sedis|complex|group|mitosporic|idae$|ineae$|ineae$|Dikarya|Eukaryota|Opisthokonta|myceta$|myceta$|Basal fungal lineages|idae$|idae$|TBM clade|^&amp|^apos|Rozella clade|mycotina$|mycotina$|Fungi|mycota$|etes$|ales$|aceae$|mitosporic|^[a-z])/) {
				#print "**$x\n";#test
					$genus = $x;
					if ($genus =~ /\w+\./) {
						$genus =~ /(\w+)\./;
						$genus = $1;
						$genus{$ref_gb} = $genus;
#					$species=();
#					$species{$ref_gb} = $species;
						$d++;
						next;
					}
					else {
						$genus{$ref_gb} = $genus;
					}
				}
			}
		}
		$d++;
		#print "genus: $genus\tspecies: $species\n";#test
	}
	$d=0;
	$e=();
	$c++;
	@line=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
#	$species=();
	$genus=();
}	

#while (my($key,$value) = each(%species)) {
#	print "key:$key\tvalue:$value\n";
#}#test

}

##########For each gi, and for each rank, compare classifiation hashes to calculate recovery and erroneous recovery
sub compare_classifications {

#declare var
my $gi;
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

#declare array

open (OUT,">>","comparison.txt") || die "Error cannot write to comparison.txt:$!\n";

while (($gb, $kingdom) = each(%kingdom) ) {	##changed $gi to $gb
	$kingdom = $kingdom{$gb};
	$ref_kingdom = $ref_kingdom{$gb};
	if ($kingdom =~ /^\S+/) {
		$coverage++;
		if ($ref_kingdom eq $kingdom) {
			$correct++;
		}
		elsif ($ref_kingdom ne $kingdom) {
			$incorrect++;
		}
	}
	$kingdom=();
	$ref_kingdom=();
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
			$correct++;
		}
		elsif ($ref_phylum ne $phylum) {
			$incorrect++;
		}
	}	
	$ref_phylum=();
	$phylum=();
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
			$correct++;
		}
		elsif ($ref_class ne $class) {
			$incorrect++;
		}
	}
	$ref_class=();
	$class=();
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
			$correct++;
		}
		elsif ($ref_order ne $order) {
			$incorrect++;
		}
	}
	$ref_order=();
	$order=();
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
			$correct++;
		}
		elsif ($ref_family ne $family) {
			$incorrect++;
		}
	}
	$ref_family=();
	$family=();
}

print OUT "family\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

while (($gb, $genus) = each(%genus) ) {	
	$ref_genus = $ref_genus{$gb};
	$genus = $genus{$gb};
	#print "ref: $ref_genus\tclassified: $genus\n";#test
	if ($genus =~ /^\S+/) {
		$coverage++;
		if ($ref_genus eq $genus) {
			$correct++;
		}
		elsif ($ref_genus ne $genus) {
			$incorrect++;
		}
	}
	$ref_genus=();
	$genus=();
}

print OUT "genus\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

#while (($gb, $species) = each(%species) ) {	
#	$ref_species = $ref_species{$gb};
#	$species = $species{$gb};
#	if ($species =~ /^\S+/) {
#		$coverage++;
#		if ($ref_species eq $species) {
#			$correct++;
#		}
#		elsif ($ref_species ne $species) {
#			$incorrect++;
#		}
#	}
#	$ref_species=();
#	$species=();
#}

#print OUT "species\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;
close OUT;

}
