#!/usr/bin/perl
#Oct. 18, 2011 edited to compare output from error_check1.plx and error_check2.plx, calc recovery and coverage wrt 0% error
#Oct.17, 2011 by Terri Porter
#Script to compare output from error_check1.plx and error_check2.plx to see difference in recovery
#usage perl error_check3.plx phylum.ref order.ref family.ref genus.ref phylum.error order.error family.error genus.error

#declare var
my $i=0;
my $line;
my $gi;
my $order;
my $ref;
my $correct=0;
my $phylum;
my $family;
my $genus;
my $coverage=0;
my $ref_coverage=0;
my $err_coverage=0;

#declare array
my @phylum_ref;
my @order_ref;
my @order_error;
my @line;
my @phylum_error;
my @family_ref;
my @family_error;
my @genus_ref;
my @genus_error;

#declare hash
my %phylum_ref;
my %order_ref;
my %family_ref;
my %genus_ref;

open (PHYLUM_REF,"<",$ARGV[0]) || die "Error cannot read from phylum.ref: $!\n";
@phylum_ref = <PHYLUM_REF>;
close PHYLUM_REF;
while ($phylum_ref[$i]) {
	$line = $phylum_ref[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi=$line[0];
	$phylum = $line[1];
	$phylum_ref{$gi} = $phylum;
	$i++;
	@line=();
	$gi=();
	$phylum=();
}
$ref_coverage = $i;
$i=0;

open (PHYLUM_ERR,"<",$ARGV[4]) || die "Error cannot read from phylum.error: $!\n";
@phylum_error = <PHYLUM_ERR>;
close PHYLUM_ERR;
while($phylum_error[$i]) {
	$line = $phylum_error[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi = $line[0];
	$phylum = $line[1];
	$ref = $phylum_ref{$gi};
	
	if ($ref) {
		if ($phylum =~ /^\S+/) {
			if ($ref eq $phylum) {
				$correct++;
				$err_coverage++;
			}
			else {
				$err_coverage++;
			}
		}
	}
	$i++;
	@line=();
	$ref=();
	$gi=();
	$phylum=();
}
$i=0;
print "Phylum correct: $correct\tPhylum coverage: $err_coverage\tRef coverage: $ref_coverage\n";
$correct=0;
$err_coverage=0;
$ref_coverage=0;


open (ORDER_REF,"<",$ARGV[1]) || die "Error cannot read from order.ref: $!\n";
@order_ref = <ORDER_REF>;
close ORDER_REF;
while ($order_ref[$i]) {
	$line = $order_ref[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi=$line[0];
	$order = $line[1];
	$order_ref{$gi} = $order;
	$i++;
	@line=();
	$gi=();
	$order=();
}
$ref_coverage = $i;
$i=0;

open (ORDER_ERR,"<",$ARGV[5]) || die "Error cannot read from order.error: $!\n";
@order_error = <ORDER_ERR>;
close ORDER_ERR;
while($order_error[$i]) {
	$line = $order_error[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi = $line[0];
	$order = $line[1];
	$ref = $order_ref{$gi};
	
	if ($ref) {
		if ($order =~ /^\S+/) {
			if ($ref eq $order) {
				$correct++;
				$err_coverage++;
			}
			else {
				$err_coverage++;
			}
		}
	}
	$i++;
	@line=();
	$ref=();
	$gi=();
	$order=();
}
$i=0;
print "Order correct: $correct\tOrder coverage: $err_coverage\tRef coverage: $ref_coverage\n";
$correct=0;
$err_coverage=0;
$ref_coverage=0;

open (FAMILY_REF,"<",$ARGV[2]) || die "Error cannot read from family.ref: $!\n";
@family_ref = <FAMILY_REF>;
close FAMILY_REF;
while ($family_ref[$i]) {
	$line = $family_ref[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi=$line[0];
	$family=$line[1];
	$family_ref{$gi} = $family;
	$i++;
	@line=();
	$gi=();
	$family=();
}
$ref_coverage=$i;
$i=0;

open (FAMILY_ERR,"<",$ARGV[6]) || die "Error cannot read from family.error: $!\n";
@family_error = <FAMILY_ERR>;
close FAMILY_ERR;
while($family_error[$i]) {
	$line = $family_error[$i];
	chomp $line;
	@line=split(/\t/,$line);
	$gi=$line[0];
	$family=$line[1];
	$ref=$family_ref{$gi};
	if ($ref) {
		if ($family =~ /^\S+/) {
			if ($ref eq $family) {
				$correct++;
				$err_coverage++;
			}
			else {
				$err_coverage++;
			}
		}
	}
	$i++;
	@line=();
	$ref=();
	$gi=();
	$family=();
}
$i=0;
print "Family correct: $correct\tFamily coverage: $err_coverage\tRef coverage: $ref_coverage\n";
$correct=0;
$err_coverage=0;
$ref_coverage=0;

open (GENUS_REF,"<",$ARGV[3]) || die "Error cannot read from genus.ref: $!\n";
@genus_ref = <GENUS_REF>;
close GENUS_REF;
while($genus_ref[$i]) {
	$line = $genus_ref[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi=$line[0];
	$genus=$line[1];
	$genus_ref{$gi}=$genus;
	$i++;
	@line=();
	$gi=();
	$genus=();
}
$ref_coverage=$i;
$i=0;

open (GENUS_ERR,"<",$ARGV[7]) || die "Error cannot read from genus.error: $!\n";
@genus_error=<GENUS_ERR>;
close GENUS_ERR;
while ($genus_error[$i]) {
	$line=$genus_error[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$gi=$line[0];
	$genus=$line[1];
	$ref=$genus_ref{$gi};
	if ($ref) {
		if ($genus =~ /^\S+/) {
			if ($ref eq $genus) {
				$correct++;
				$err_coverage++;
			}
			else {
				$err_coverage++;
			}
		}
	}
	$i++;
	@line=();
	$ref=();
	$gi=();
	$genus=();
}
$i=0;
print "Genus correct: $correct\tGenus coverage: $err_coverage\tRef coverage: $ref_coverage\n";
$correct=0;
$err_coverage=0;
$ref_coverage=0;
	
