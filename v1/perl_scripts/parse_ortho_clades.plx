#!/usr/bin/perl
#Dec. 16/11 by Terri Porter
#Script to parse ortho.clades from parse_phylip_outfiles.plx
#usage perl parse_ortho_clades.plx ortho.clades

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $id_current;
my $id_next;
my $pattern;
my $support;
my $tot_boot;
my $percent;
my $count_support;
my $return_support=0;
my $j=0;
my $flag=0;
my $k=0;
my $next_line;
my $basid_support;
my $tot_orthos;
my $single_taxon_clade;
my $return_support2;
my $basid_single_taxon_clade;

#declare array
my @in;
my @line;
my @pattern;
my @reference_taxa = ("LB","UM","SR","AN","SC","AA","MV","EP","RO","CR","NE","SP","BD","PD","CA","BB","BE","AM","MB");
my @count_support;
my @single_taxon_clade;
my @next_line;

#declare hash
my %pattern;
my %basid_support;
my %basid_single_taxon_clade;
my %orthos;

open (IN,"<", $ARGV[0]) || die "Error cannot read ortho.clades: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	$j=$i+1;

	$next_line = $in[$j];
	
	if ($next_line) {
		chomp $next_line;
		@next_line = split(/\t/,$next_line);
		$id_next = $next_line[0];
		#print "got next_line\n";
	}
	elsif (!$next_line) {
		#print "found eof\n";
		$next_line = "eof";
		$id_next = "nil";
	}
	
	if ($i==0) {
		$i++;
		next;
	}
	else {
		@line = split(/\t/,$line);
		$id_current = $line[0];
		$orthos{$id_current}=1;
		$pattern = $line[1];
		$support = $line[2];
		$tot_boot = $line[3];
		$percent = $support/$tot_boot;
		#print "percent: $percent\n";
		@pattern = split(//,$pattern);
		$pattern{LB} = $pattern[0];
		$pattern{UM} = $pattern[1];
		$pattern{SR} = $pattern[2];
		$pattern{AN} = $pattern[3];
		$pattern{SC} = $pattern[4];
		$pattern{AA} = $pattern[5];
		$pattern{MV} = $pattern[6];
		$pattern{EP} = $pattern[7];
		$pattern{RO} = $pattern[8];
		$pattern{CR} = $pattern[9];
		$pattern{NE} = $pattern[10];
		$pattern{SP} = $pattern[11];
		$pattern{BD} = $pattern[12];
		$pattern{PD} = $pattern[13];
		$pattern{CA} = $pattern[14];
		$pattern{BB} = $pattern[15];
		$pattern{BE} = $pattern[16];
		$pattern{AM} = $pattern[17];
		$pattern{MB} = $pattern[18];

		#print "id_current: $id_current\tid_next: $id_next\n";
		#test
#		while ( my ($key,$value) = each (%pattern) ) {
#			print "$key => $value\n";
#		}
		if ($percent >= 0.700) {
			#check clade support
			#print "got > 70% clade support\n";
			check_basid();
			assess_support();
			if ($return_support == 1) {
				$basid_support{$id_current} = 1;
			}
			$return_support=0;
			if ($return_support ==2) {
				$basid_single_taxon_clade{$id_current} = 1;
			}
	#			check_asco();
#			check_mucoro();
#			check_blastoclad();
#			check_chytrid();
			#check relationship support
#			check_dikarya();
#			check_dikarya_mucoro();
#			check_dikarya_mucoro_blastoclad();
#			check_dikarya_mucoro_blastoclad_chytrid();
			##check alternative support
#			check_dikarya_mucoro_chytrid();
		}
	}	
	$i++;
	$id_current=();
	$id_next=();
	$pattern=();
	$support=();
	$tot_boot=();
	$percent=();
	@line=();
	@pattern=();
	%pattern=();
	@next_line=();
}
$i=0;

$basid_support = keys(%basid_support);
print "basid_support = $basid_support\n";
$basid_single_taxon_clade = keys(%basid_single_taxon_clade);
print "basid_single_taxon_clade = $basid_single_taxon_clade\n";

$tot_orthos = keys(%orthos);
print "total number of orthos: $tot_orthos\n";

####################

sub assess_support {

	if ($id_current eq $id_next) {
		#print "id_current eq id_next\n";
		push(@count_support,$count_support);
		push(@single_taxon_clade,$single_taxon_clade);
	}
	elsif ($id_current ne $id_next) {
		#print "id_current ne id_next\n";
		while ($count_support[$k]) {
			$count_support = $count_support[$k];
			if ($flag==0 && $count_support >=1) {
				$flag=1;
			}
			elsif ($flag==0 && $count_support<1) {
				$k++;
				next;
			}
			elsif ($flag==1) {
				$k++;
				next;
			}
			$k++;
		}
		$k=0;

		if ($flag==1) {
			$return_support = 1;
		}
		$flag=0;

		while ($single_taxon_clade[$k]) {
			$single_taxon_clade = $single_taxon_clade[$k];
			if ($flag==0 && $single_taxon_clade >=1) {
				$flag=1;
			}
			elsif ($flag==0 && $single_taxon_clade<1) {
				$k++;
				next;
			}
			elsif ($flag==1) {
				$k++;
				next;
			}
			$k++;
		}
		$k=0;

		if ($flag==1) {
			$return_support2 = 1;
		}
		$flag=0;

		@count_support=();
		@single_taxon_clade=();
		#print $count_support."\n";
		push(@count_support,$count_support);
		push(@single_taxon_clade,$single_taxon_clade);
	}
}

####################

sub check_basid {

	$count_support=0;
	$single_taxon_clade=0;

	if ($pattern{LB} eq "*" && $pattern{UM} eq "*" && $pattern{SR} eq "*") {
		$count_support++;
	}
	elsif ($pattern{LB} eq "*" && $pattern{UM} eq "*" && $pattern{SR} eq "-") {
		$count_support++;
	}
	elsif ($pattern{LB} eq "*" && $pattern{UM} eq "-" && $pattern{SR} eq "*") {
		$count_support++;
	}
	elsif ($pattern{LB} eq "-" && $pattern{UM} eq "*" && $pattern{SR} eq "*") {
		$count_support++;
	}
	elsif ($pattern{LB} eq "*" && $pattern{UM} eq "-" && $pattern{SR} eq "-") {
		$single_taxon_clade++;
	}
	elsif ($pattern{LB} eq "-" && $pattern{UM} eq "*" && $pattern{SR} eq "-") {
		$single_taxon_clade++;
	}
	elsif ($pattern{LB} eq "-" && $pattern{UM} eq "-" && $pattern{SR} eq "*") {
		$single_taxon_clade++;
	}


}

####################
