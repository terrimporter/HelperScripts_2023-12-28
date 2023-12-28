#!/usr/bin/perl
#March 14, 2012 edited to print gi\tname\n format for error checking, copied get_recovery3b.plx and get_recovery3_error1.plx
#NEW USAGE perl get_recovery3b_error1.plx features.txt.parsed gb.query gi.query NBC_download.txt

#Sept.30, 2011 edited to parse NBC output and filter results by statistical cutoff (80 default)
#new usage perl get_recovery2.plx features.txt.parsed NBC_download.fasta
#Sept.29, 2011 by Terri Porter
#Script to compare features.txt.parsed from test set with megan_hit_lineage.txt.mapped
#usage perl get_recovery.plx features.txt.parsed gb.query gi.query megan_hit_lineage.txt.mapped

use strict;
use warnings;

#declare var
my $CUTOFF=50;###change to 50 for <250bp, otherwise 80 cutoff###

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

open (GB,"<",$ARGV[1]) || die "Error cannot open gb.query:$!\n";
@gb=<GB>;
close GB;

open (GI,"<",$ARGV[2]) || die "Error cannot open gi.query:$!\n";
@gi=<GI>;
close GI;

open (NBC,"<",$ARGV[3]) || die "Error cannot open nbc.fasta_download.txt:$!\n";
@nbc = <NBC>;
close NBC;

make_map();

parse_ref_lineage();

parse_NBC_lineage();

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
#	print "map key:$key\tvalue:$value\n";
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

#while(my($key,$value) = each(%ref_genus)) {
#	print "ref_genus key:$key\tvalue:$value\n";
#}#test

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
	if ($line =~ /^Details:/) {
		$flag=1;
	}
	if ($flag ==1) {
		@line = split(/;/,$line);
		$ref_gb = shift @line;
		#@fragment1 = split(/\s+/,$fragment1);
		#$ref_gi = $fragment1[0];
	
		while ($line[$d]) { #iterate through lineage items (assignments and probabilities)
			$x = $line[$d];
			$x =~ s/^\s{1}//;
			if ($x =~ /&amp;apos;/) {
				$x =~ s/&amp;apos;//;
			}

			if ($x eq "Fungi") {
				$f="nil";
				$kingdom = $x;
				$kingdom{$ref_gb} = $kingdom;
				$e = $d+1;
				$stat = $line[$e];
				$stat =~ /^\s{1}(\d+)%/;
				$prob = $1;
				$kingdom_stat{$ref_gb} = $prob;
				$e=();
				$stat=();
				$prob=();
			}
			elsif ($x =~ /mycota$/) {
				$f="nil";
				if ($x !~ /(mitosporic|Euglenida)/) {
					$phylum =$x;
					$phylum{$ref_gb} = $phylum;
					$e = $d+1;
					$stat = $line[$e];
					$stat =~ /^\s{1}(\d+)%/;
					$prob = $1;
					$phylum_stat{$ref_gb} = $prob;
					$e=();
					$stat=();
					$prob=();
				}
			}
			elsif ($x =~ /etes$/) {
				$f="nil";
				$class =$x;
				$class{$ref_gb} = $class;
				$e = $d+1;
				$stat = $line[$e];
				$stat =~ /^\s{1}(\d+)%/;
				$prob = $1;
				$class_stat{$ref_gb} = $prob;
				$e=();
				$stat=();
				$prob=();
			}
			elsif ($x =~ /ales$/) {
				$f="nil";
				$order =$x;
				if ($x !~ /(mitosporic|Choanoflagellida)/) {
					$order{$ref_gb} = $order;
					$e = $d+1;
					$stat = $line[$e];
					$stat =~ /^\s{1}(\d+)%/;
					$prob = $1;
					$order_stat{$ref_gb} = $prob;
					$e=();
					$stat=();
					$prob=();
				}
			}
			elsif ($x =~ /aceae$/) {
				$f="nil";
				$family =$x;
				if ($x !~ /(mitosporic|Codonosigidae)/) {
					$family{$ref_gb} = $family;
					$e = $d+1;
					$stat = $line[$e];
					$stat =~ /^\s{1}(\d+)%/;
					$prob = $1;
					$family_stat{$ref_gb} = $prob;
					$stat=();
					$prob=();
#					$f=$e+1;#genus index  ####this is a BAAAD way to find the genus!!!
#					$g=$e+2;#genus_stat index
					$e=();
				}
			}
#			elsif ($f ne "nil") {
#				$genus = $line[$f];
#				$genus =~ s/^\s{1}//;
#				if ($genus !~ /(incertae sedis|complex|group)/) {
#					$genus{$ref_gb} = $genus;
#					if ($g ne "nil") {
#						$stat = $line[$g];
#						$stat =~ /^\s+(\d+)%/;
#						$prob = $1;
#						$genus_stat{$ref_gb} = $prob;
#						$e=();
#						$f="nil";
#						$g="nil";
#						$stat=();
#						$prob=();
#					}
#				}
#			}
			$d++;
		}
		$d=0;
		$e=();
#		$f="nil";
#		$g="nil";
#		@line=();
		$kingdom=();
		$phylum=();
		$class=();
		$order=();
		$family=();
#		$species=();
#		$genus=();
		
#		print "array: @line\n";#test
		$g = scalar(@line); #genus stat index
#		print "g_index: $g\t";
		$g=$g-1;		
		my $test1 = $line[$g];
#		print "g_value: $test\n";
		
		$f=$g-1; #genus index
#		print "f_index: $f\t";
		
		my $test2 = $line[$f];
#		print "f_value: $test\n";
	
		if ($line[$f]) {
			$genus = $line[$f];
			$genus =~ s/^\s{1}//;
			if ($genus !~ /(incertae sedis|complex|group|mitosporic|Fungi|mycota$|etes$|ales$|aceae$)/) {
				$genus{$ref_gb} = $genus;
#				print "$genus\t$test1\n";
				if ($g ne "nil") {
					$stat = $line[$g];
					$stat =~ /^\s+(\d+)%/;
					$prob = $1;
					$genus_stat{$ref_gb} = $prob;
#					$e=();
#					$f="nil";
#					$g="nil";
#					$stat=();
#					$prob=();
				}
			}
		}
		$f="nil";
		$g="nil";
		$test1=();
		$test2=();
		$stat=();
		$prob=();
		@line=();

	}
	$c++;
}	
#while(my($key,$value) = each(%genus)) {
#	print "genus key:$key\tvalue:$value\n";
#}#test

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
my $gi="nil";

#declare array

#open (OUT,">>","comparison.txt") || die "Error cannot write to comparison.txt:$!\n";

#open (CORRECT,">>","correct.txt") || die "Error cannot write to correct.txt: $!\n";

#open (INCORRECT,">>","incorrect.txt") || die "Error cannot write to incorrect.txt: $!\n";

while (($gb, $kingdom) = each(%kingdom) ) {	
	$kingdom = $kingdom{$gb};
	$ref_kingdom = $ref_kingdom{$gb};
	if ($kingdom =~ /^\S+/) {
		$coverage++;
		if ($ref_kingdom eq $kingdom) {
			$kingdom_stat = $kingdom_stat{$gb};
			if ($kingdom_stat >= $CUTOFF) {
				$correct++;
#				print CORRECT "$gb\t$kingdom\t$kingdom_stat\n";
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_kingdom ne $kingdom) {
			if ($kingdom_stat >= $CUTOFF) {
				$incorrect++;
#				print INCORRECT "$gb\t$kingdom\t$kingdom_stat\n";
			}
			else {
				$coverage--;
			}
		}
	}
}

#print OUT "rank\tcorrect\tincorrect\tcoverage\n";
#print OUT "kingdom\t$correct\t$incorrect\t$coverage\n";
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
#				print CORRECT "$gb\t$phylum\t$phylum_stat\n";
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_phylum ne $phylum) {
			if ($phylum_stat >= $CUTOFF) {
				$incorrect++;
#				print INCORRECT "$gb\t$phylum\t$phylum_stat\n";
			}
			else {
				$coverage--;
			}
		}
	}
}

#print OUT "phylum\t$correct\t$incorrect\t$coverage\n";
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
#				print CORRECT "$gb\t$class\t$class_stat\n";
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_class ne $class) {
			if ($class_stat >= $CUTOFF) {
				$incorrect++;
#				print INCORRECT "$gb\t$class\t$class_stat\n";
			}
			else {
				$coverage--;
			}
		}
	}
}

#print OUT "class\t$correct\t$incorrect\t$coverage\n";
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
#				print CORRECT "$gb\t$order\t$order_stat\n";
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_order ne $order) {
			if ($order_stat >= $CUTOFF) {
				$incorrect++;
#				print INCORRECT "$gb\t$order\t$order_stat\n";
			}
			else {
				$coverage--;
			}
		}
	}
}

#print OUT "order\t$correct\t$incorrect\t$coverage\n";
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
#				print CORRECT "$gb\t$family\t$family_stat\n";
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_family ne $family) {
			if ($family_stat >= $CUTOFF) {
				$incorrect++;
#				print INCORRECT "$gb\t$family\t$family_stat\n";
			}
			else {
				$coverage--;
			}
		}
	}
}

#print OUT "family\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

#edit genus rank only to print gi\tname\n formatted table

open (TABLE,">>", "genus.ref") || die "Error cannot write to genus.ref: $!\n";

while (($gb, $genus) = each(%genus) ) {	
	$ref_genus = $ref_genus{$gb};
	$genus = $genus{$gb};
	$gi = $map_gb_gi{$gb};
	#print "Ref genus: $ref_genus\tgenus: $genus\n";#test
	if ($genus =~ /^\S+/) {
		$coverage++;
		if ($ref_genus eq $genus) {
			$genus_stat = $genus_stat{$gb};
			if ($genus_stat >= $CUTOFF) {
				$correct++;
#				print "$gb\t$genus\t$genus_stat\n";
				print TABLE "$gi\t$genus\n";
			}
			else {
				$coverage--;
			}
		}
		elsif ($ref_genus ne $genus) {
			if ($genus_stat >= $CUTOFF) {
				$incorrect++;
#				print  INCORRECT "$gb\t$genus\t$genus_stat\n";
			}
			else {
				$coverage--;
			}
		}
	}
	$gi="nil";
}

#print OUT "genus\t$correct\t$incorrect\t$coverage\n";
$coverage=0;
$correct=0;
$incorrect=0;

#while (($gb, $species) = each(%species) ) {	
#	$ref_species = $ref_species{$gb};
#	$species = $species{$gb};
#	if ($species !~ /^\s+/) {
#		$coverage++;
#	}
#	if ($ref_species eq $species) {
#		$correct++;
#	}
#	elsif ($ref_species ne $species) {
#		$incorrect++;
#	}
#}

#print OUT "species\t$correct\t$incorrect\t$coverage\n";
#$coverage=0;
#$correct=0;
#$incorrect=0;
#close OUT;

close TABLE;

}
