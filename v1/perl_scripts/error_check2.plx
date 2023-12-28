#!/usr/bin/perl
#Oct. 18, 2011 edited to process simulated error results into gi\trank\n format
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
	$ref_gi = $map_gb_gi{$ref_gb};
	
	$ref_organism = $line[1];
	$ref_organism{$ref_gi}= $ref_organism;
	@organism = split(/ /,$ref_organism);
	$ref_genus = $organism[0];
	$ref_genus{$ref_gi} = $ref_genus;
	$ref_species = $organism[1];
	$ref_species{$ref_gi} = $ref_species;
	
	$ref_strain = $line[2];
	$ref_isolate = $line[3];
	$ref_length = $line[4];
	$ref_classification = $line[5];
	@classification = split(/ /,$ref_classification);

	foreach $x (@classification) {
		if ($x eq "Fungi") {
			$ref_kingdom = $x;
			$ref_kingdom{$ref_gi} = $ref_kingdom;
		}
		elsif ($x =~ /mycota$/) {
			$ref_phylum = $x;
			$ref_phylum{$ref_gi} = $ref_phylum;
		}
		elsif ($x =~ /etes$/) {
			$ref_class = $x;
			$ref_class{$ref_gi} = $ref_class;
		}
		elsif ($x =~ /ales$/) {
			$ref_order = $x;
			$ref_order{$ref_gi} = $ref_order;
		}
		elsif ($x =~ /aceae$/) {
			$ref_family = $x;
			$ref_family{$ref_gi} = $ref_family;
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
	$ref_gi = $fragment1[0];
	
	while ($line[$d]) {
		$x = $line[$d];
		$x =~ s/^\s{1}//;
		$x =~ s/^&amp;//g;
		$x =~ s/^apos;//g;
		
		if ($x =~ "Fungi") {
			if ($x =~ /Fungi$/) {
				$kingdom = $x;
				$kingdom{$ref_gi} = $kingdom;
			}
			elsif ($x =~ /Fungi\./) {
				$x =~ /(Fungi)\.$/;
				$kingdom = $1;
				$kingdom{$ref_gi} = $kingdom;
				$d++;
				next;
			}
			#print "kingdom: $kingdom\n";#test
		}
		elsif ($x =~ /mycota/) {
			if ($x !~ /(mitosporic|Euglenida)/) {
				if ($x =~ /mycota$/) {
					$phylum =$x;
					$phylum{$ref_gi} = $phylum;
				}
				elsif ($x =~ /mycota\.$/) {
					$x =~ s/\.$//;
					$phylum = $x;
					$phylum{$ref_gi} =$phylum;
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
					$class{$ref_gi} = $class;
				}
				elsif ($x =~ /etes\.$/) {
					$x =~ s/\.$//;
					$class =$x;
					$class{$ref_gi} = $class;
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
					$order{$ref_gi} = $order;
				}
				elsif ($x =~ /(ales)\.$/) {
					$x =~ s/\.$//;
					$order =$x;
					$order{$ref_gi} = $order;
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
					$family{$ref_gi} = $family;
				}
				elsif ($x =~ /aceae\.$/) {
					$x =~ s/\.$//;
					$family =$x;
					$family{$ref_gi} = $family;
					$d++;
					$genus =();
					$genus{$ref_gi} = $genus;
					$species = ();
					$species{$ref_gi} = $species;
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
		elsif ($x =~ /\.\.$/) {
			$x =~ s/^\.+//;
			$x =~s/^\s+//;
			$x =~ s/\.+$//;
			$species = $x;
			$species{$ref_gi} = $species;
			#print "species:$species\n";
		}
		else {
			if ($x !~ /(incertae sedis|complex|group|mitosporic|idae$|ineae$|ineae\.$|Dikarya|Eukaryota|Opisthokonta|myceta$|myceta\.$|Basal fungal lineages|idae$|idae\.$|TBM clade|^&amp|^apos|Rozella clade|mycotina$|mycotina\.$)/) {
				#print "**$x\n";#test
				$genus = $x;
				if ($genus =~ /\w+\./) {
					$genus =~ /(\w+)\./;
					$genus = $1;
					$genus{$ref_gi} = $genus;
					$species=();
					$species{$ref_gi} = $species;
					$d++;
					next;
				}
				else {
					$genus{$ref_gi} = $genus;
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
	$species=();
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

open (PHYLUM,">>","phylum.error") || die "Error cannot write to phylum.error: $!\n";

while (($gi, $phylum) = each(%phylum) ) {	
	$ref_phylum = $ref_phylum{$gi};
	$phylum = $phylum{$gi};
	if ($phylum =~ /^\S+/) {
		$coverage++;
		if ($ref_phylum eq $phylum) {
			$correct++;
			print PHYLUM "$gi\t$phylum\n";
		}
		elsif ($ref_phylum ne $phylum) {
			$incorrect++;
			#print ORDER "$gi\t$order\n";
		}
	}
	$ref_phylum=();
	$phylum=();
}
close PHYLUM;
print "Phylum correct: $correct\t Phylum coverage: $coverage\n";

$coverage=0;
$correct=0;
$incorrect=0;


open (ORDER,">>","order.error") || die "Error cannot write to order.error: $!\n";

while (($gi, $order) = each(%order) ) {	
	$ref_order = $ref_order{$gi};
	$order = $order{$gi};
	if ($order =~ /^\S+/) {
		$coverage++;
		if ($ref_order eq $order) {
			$correct++;
			print ORDER "$gi\t$order\n";
		}
		elsif ($ref_order ne $order) {
			$incorrect++;
			#print ORDER "$gi\t$order\n";
		}
	}
	$ref_order=();
	$order=();
}
close ORDER;
print "Order correct: $correct\t Order coverage: $coverage\n";

$coverage=0;
$correct=0;
$incorrect=0;

open (FAMILY,">>","family.error") || die "Error cannot write to family.error: $!\n";

while (($gi, $family) = each(%family) ) {	
	$ref_family = $ref_family{$gi};
	$family = $family{$gi};
	if ($family =~ /^\S+/) {
		$coverage++;
		if ($ref_family eq $family) {
			$correct++;
			print FAMILY "$gi\t$family\n";
		}
		elsif ($ref_family ne $family) {
			$incorrect++;
			print FAMILY "$gi\t$family\n";
		}
	}
	$ref_family=();
	$family=();
}
close FAMILY;
print "Family correct: $correct\t Family coverage: $coverage\n";

$coverage=0;
$correct=0;
$incorrect=0;

open (GENUS,">>","genus.error") || die "Error cannot write to genus.error: $!\n";

while (($gi, $genus) = each(%genus) ) {	
	$ref_genus = $ref_genus{$gi};
	$genus = $genus{$gi};
	if ($genus =~ /^\S+/) {
		$coverage++;
		if ($ref_genus eq $genus) {
			$correct++;
			print GENUS "$gi\t$genus\n";
		}
		elsif ($ref_genus ne $genus) {
			$incorrect++;
			print GENUS "$gi\t$genus\n";
		}
	}
	$ref_genus=();
	$genus=();
}
close GENUS;

print "Genus correct: $correct\t Genus coverage: $coverage\n";

$coverage=0;
$correct=0;
$incorrect=0;

close OUT;

}

####################


