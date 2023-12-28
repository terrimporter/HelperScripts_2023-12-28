#!/usr/bin/perl
# Teresita M. Porter, Aug. 20, 2020
# Script to process SILVA records to genus rank only
# USAGE perl process_SILVA_138.plx tax_slv_ssu_138.txt SILVA_138_SSURef_NR99_tax_silva.fasta

use strict;
use warnings;
use Data::Dumper;

# vars
my $i=0;
my $j;
my $line;
my $taxonline;
my $taxonpath;
my $prefix_old;
my $taxon;
my $taxid;
my $rank;
my $header;
my $seq;
my $flag=0; # flag
my $first; # first part of the fasta header with '>' and accession
my $newtaxonpath;
my $second; # second part of fasta header with taxonpath
my $prefix;
my $counter=1;
my $prevtermcounter;
my $taxline;
my $domain;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;

# arrays
my @tax;
my @line;
my @taxonline;
my @taxonpath;
my @fas;
my @header;
my @second; # taxonpath in parts
my @taxline;
my @taxon;

# hashes
my %globaltax; # key = taxonpath, value = rank
my %recordtax; # key = taxon, value=rank
my %tax; # key = taxon, value=counter*taxon*prevtaxoncounter*rank#*rank
my %uniques; # key = taxon, value = lineage prefix up to but not including current taxon based on SILVA prefixes
my %duplicates; # key = taxon, value = 1 based on SILVA prefixes
#my %uniques2; # as above but just including whole ranks
#my %duplicates2;

# hash taxonomy for easy rank lookup
open (IN, "<", $ARGV[0]) || die "Cannot open taxonomy file:$!\n";
@tax = <IN>;
close IN;

while($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$taxonpath = $line[0];
	$taxonpath =~ s/ /_/g;
	$taxonpath =~ s/-/_/g;

	$taxid = $line[1];

	$rank = $line[2];

	$globaltax{$taxonpath} = $rank;

	# keep track of duplicate taxa
	@taxonpath = split(';', $taxonpath);

	$j = scalar(@taxonpath);

	while ($j > 0) {
		$taxon = pop(@taxonpath);
		$prefix = join(';',@taxonpath);

		if (!exists $uniques{$taxon}) {
			$uniques{$taxon} = $prefix;
		}
		else {
			$prefix_old = $uniques{$taxon};
			
			if ($prefix_old ne $prefix) {
				$duplicates{$taxon} = 1;
			}
		}
		$j--;
	}

	$i++;

}
$i=0;

# for each record, go through fill in missing ranks and check for dups
open (IN, "<", $ARGV[1]) || die "Cannot open fasta file: $!\n";
@fas = <IN>;
close IN;

while($fas[$i]) {
	$line = $fas[$i];
	chomp $line;

	# start of new record
	if ($line =~ /^>/) {

		# process previous record with reformatted header
		if (length($seq)) {
			if ($recordtax{'genus'}) {
				# print to TAX file
#				process_tax();
				%recordtax=();
				$seq=();

			}
		}

		# process current header
		$flag = 0;
		$header = $line;
		@header = split(' ',$header);
		$first = shift(@header);
		$second = join('_',@header);

		# build up a record hash
		build_hash($second);

		# ensure we can find the domain, else go on to next record
		if (!exists $recordtax{"domain"}){
			$i++;
			%recordtax=();
			next;
		}

		# ensure we can find the kingdom 
		if (!exists $recordtax{"kingdom"}) {
			$prefix = $recordtax{"domain"};
			$taxon = $prefix."_undef";
			$recordtax{"kingdom"} = $taxon;
		}

		# ensure we can find the phylum, check if new 'undef' taxon is a duplicate
		if (!exists $recordtax{"phylum"}) {
			$prefix = $recordtax{"kingdom"};
			$taxon = $prefix."_undef";
			$recordtax{"phylum"} = $taxon;

			if (!exists $uniques{$taxon}) {
				$prefix = $recordtax{"domain"}.";".$recordtax{"kingdom"};
				$uniques{$taxon} = $prefix;
			}
			else {
				$prefix_old = $uniques{$taxon};
				
				if ($prefix_old ne $prefix) {
					$duplicates{$taxon} = 1;
				}
			}
		}


		# ensure we can find the class, check if new 'undef' taxon is a duplicate
		if (!exists $recordtax{"class"}) {
			$prefix = $recordtax{"phylum"};
			$taxon = $prefix."_undef";
			$recordtax{"class"} = $taxon;
			if (!exists $uniques{$taxon}) {
				$prefix = $recordtax{"domain"}.";".$recordtax{"kingdom"}.";".$recordtax{"phylum"};
				$uniques{$taxon} = $prefix;
			}
			else {
				$prefix_old = $uniques{$taxon};
				
				if ($prefix_old ne $prefix) {
					$duplicates{$taxon} = 1;
				}
			}
		}
	

		# ensure we can find the order, check if new 'undef' taxon is a duplicate
		if (!exists $recordtax{"order"}) {
			$prefix = $recordtax{"class"};
			$taxon = $prefix."_undef";
			$recordtax{"order"} = $taxon;
		
			if (!exists $uniques{$taxon}) {
				$prefix = $recordtax{"domain"}.";".$recordtax{"kingdom"}.";".$recordtax{"phylum"}.";".$recordtax{"class"};
				$uniques{$taxon} = $prefix;
			}
			else {
				$prefix_old = $uniques{$taxon};
			
				if ($prefix_old ne $prefix) {
					$duplicates{$taxon} = 1;
				}
			}
		}
	

		# ensure we can find the family, check if new 'undef' taxon is a duplicate
		if (!exists $recordtax{"family"}) {
			$prefix = $recordtax{"order"};
			$taxon = $prefix."_undef";
			$recordtax{"family"} = $taxon;
	
			if (!exists $uniques{$taxon}) {
				$prefix = $recordtax{"domain"}.";".$recordtax{"kingdom"}.";".$recordtax{"phylum"}.";".$recordtax{"class"}.";".$recordtax{"order"};
				$uniques{$taxon} = $prefix;
			}
			else {
				$prefix_old = $uniques{$taxon};
			
				if ($prefix_old ne $prefix) {
					$duplicates{$taxon} = 1;
				}
			}
		}
	

		# ensure we can find the genus, if not go on to next record
		if (!exists $recordtax{"genus"}) {
			%recordtax=();
			$i++;
			next;
		}

	}
	# start of new seq
	elsif ($flag == 0) {
		$flag = 1;
		$seq = $line;
	}
	# continuation of seq on next line
	elsif ($flag == 1) {
		$seq = $seq.$line;
	}
	$i++;

}
$i=0;

# process last record
if (length($seq)) {
	if ($recordtax{'genus'}) {

		# print to TAX file
#		process_tax();
		%recordtax=();
		$seq=();
	}
}

# reformat fasta for RDP, drop species, ensure an entry for each rank even if missing from SILVA
open (FAS, ">>", "testNBC.fasta") || die "Cannot open fasta output file: $!\n";
open (TAX, ">>", "testNBC.taxonomy") || die "Cannot open taxonomy output file: $!\n";

while($fas[$i]) {
	$line = $fas[$i];
	chomp $line;

	# start of new record
	if ($line =~ /^>/) {

		# print out previous record with reformatted header
		if (length($seq)) {
			if ($recordtax{'genus'}) {

				# print to TAX file
				print_to_tax();
	
				# print to FASTA file
				$second = "cellularOrganisms;$recordtax{'domain'};$recordtax{'kingdom'};$recordtax{'phylum'};$recordtax{'class'};$recordtax{'order'};$recordtax{'family'};$recordtax{'genus'}";
				$seq =~ s/U/T/g;
				print FAS "$first\t$second\n$seq\n";

				%recordtax=();
				$seq=();

			}
		}

		# process current header
		$flag = 0;
		$header = $line;
		@header = split(' ',$header);
		$first = shift(@header);
		$second = join('_',@header);

		# build up a record hash
		build_hash($second);

		# ensure we can find the domain, else go on to next record
		if (!exists $recordtax{"domain"}){
			$i++;
			%recordtax=();
			next;
		}

		# ensure we can find the kingdom 
		if (!exists $recordtax{"kingdom"}) {
			$prefix = $recordtax{"domain"};
			$taxon = $prefix."_undef";
			$recordtax{"kingdom"} = $taxon;
		}

		# ensure we can find the phylum
		if (!exists $recordtax{"phylum"}) {
			$prefix = $recordtax{"kingdom"};
			$taxon = $prefix."_undef";
			$recordtax{"phylum"} = $taxon;
		}


		# ensure we can find the class
		if (!exists $recordtax{"class"}) {
			$prefix = $recordtax{"phylum"};
			$taxon = $prefix."_undef";
			$recordtax{"class"} = $taxon;
		}
	

		# ensure we can find the order
		if (!exists $recordtax{"order"}) {
			$prefix = $recordtax{"class"};
			$taxon = $prefix."_undef";
			$recordtax{"order"} = $taxon;
		}
	

		# ensure we can find the family
		if (!exists $recordtax{"family"}) {
			$prefix = $recordtax{"order"};
			$taxon = $prefix."_undef";
			$recordtax{"family"} = $taxon;
		}
	

		# ensure we can find the genus, if not go on to next record
		if (!exists $recordtax{"genus"}) {
			%recordtax=();
			$i++;
			next;
		}
	}
	# start of new seq
	elsif ($flag == 0) {
		$flag = 1;
		$seq = $line;
	}
	# continuation of seq on next line
	elsif ($flag == 1) {
		$seq = $seq.$line;
	}
	$i++;

}
$i=0;

# process last record
if (length($seq)) {
	if ($recordtax{'genus'}) {

		# print to TAX file
		print_to_tax();

		# print to FASTA file
		$second = "cellularOrganisms;$recordtax{'domain'};$recordtax{'kingdom'};$recordtax{'phylum'};$recordtax{'class'};$recordtax{'order'};$recordtax{'family'};$recordtax{'genus'}";
		$seq =~ s/U/T/g;
		print FAS "$first $second\n$seq\n";

		%recordtax=();
	}
}

close FAS;
close TAX;

print Dumper(\%duplicates);
#print Dumper(\%duplicates2);


########################################
sub build_hash{

	$second = $_[0];

	while(length($second)) {
		$second = $second.";";
		$second =~ s/ /_/g;
		$second =~ s/-/_/g;

		if (exists $globaltax{$second}) {
			$rank = $globaltax{$second};
		}
		else {
			$rank = "species";
		}
		@second = split(';',$second);
		$taxon = pop(@second);
		$recordtax{$rank} = $taxon;
		$second = join(';',@second);

	}
}

########################################
sub print_to_tax{
	if (!exists $tax{"cellularOrganisms"}){
		$tax{"cellularOrganisms"} = "$counter*cellularOrganisms*0*0*root";
		$prevtermcounter=$counter;
		$counter++;
		print TAX $tax{"cellularOrganisms"}."\n";
	}
	else {
		$taxline = $tax{"cellularOrganisms"};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$domain = $recordtax{'domain'};
	if (!exists $tax{$domain}) {
		$tax{$domain} = "$counter*$domain*$prevtermcounter*1*domain";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$domain};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$domain};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$kingdom = $recordtax{'kingdom'};
	if (!exists $tax{$kingdom}) {
		$tax{$kingdom} = "$counter*$kingdom*$prevtermcounter*2*kingdom";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$kingdom};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$kingdom};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$phylum = $recordtax{'phylum'};
	if (exists $duplicates{$phylum}) {
		$recordtax{'phylum'} = 	$kingdom."_".$phylum;
		$phylum = $recordtax{'phylum'};
	}
	if (!exists $tax{$phylum}) {
		$tax{$phylum} = "$counter*$phylum*$prevtermcounter*3*phylum";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$phylum};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$phylum};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$class = $recordtax{'class'};
	if (exists $duplicates{$class}) {
		$recordtax{'class'} = 	$phylum."_".$class;
		$class = $recordtax{'class'};
	}
	if (!exists $tax{$class}) {
		$tax{$class} = "$counter*$class*$prevtermcounter*4*class";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$class};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$class};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$order = $recordtax{'order'};
	if (exists $duplicates{$order}) {
		$recordtax{'order'} = 	$class."_".$order;
		$order = $recordtax{'order'};
	}
	if (!exists $tax{$order}) {
		$tax{$order} = "$counter*$order*$prevtermcounter*5*order";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$order};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$order};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$family = $recordtax{'family'};
	if (exists $duplicates{$family}) {
		$recordtax{'family'} = 	$order."_".$family;
		$family = $recordtax{'family'};
	}
	if (!exists $tax{$family}) {
		$tax{$family} = "$counter*$family*$prevtermcounter*6*family";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$family};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$family};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

	$genus = $recordtax{'genus'};
	if (exists $duplicates{$genus}) {
		$recordtax{'genus'} = $family."_".$genus;
		$genus = $recordtax{'genus'};
	}
	elsif ($genus eq "marine_group") { # for some reason RDP classifier doesn't parse the alphanumeric groups before marine_group properly?
		$recordtax{'genus'} = $family."_".$genus;
		$genus = $recordtax{'genus'};
	}
	if (!exists $tax{$genus}) {
		$tax{$genus} = "$counter*$genus*$prevtermcounter*7*genus";
		$prevtermcounter=$counter;
		$counter++;
		$taxline = $tax{$genus};
		print TAX $taxline."\n";
	}
	else {
		$taxline = $tax{$genus};
		@taxline = split(/\*/,$taxline);
		$prevtermcounter = shift(@taxline);
	}

}
