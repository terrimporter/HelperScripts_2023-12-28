#!/usr/bin/perl
#Sept.28, 2010 by Terri Porter
#Script to automatically generate infiles for SSA from a completely formatted SAP_whole.txt file and a map.txt that contains all properly formatted taxa.
#usage $perl make_SSA_infiles.plx map.txt SAP_whole.txt

use strict;
use warnings;

#declare variables
my $line;
my $num_taxa;
my $num_char;
my $num_species;

#declare arrays
my @line;
my @species;#species ids from map file
my @taxa_ids;
my @whole;
my @species_ids;#species ids from whole file
my @species_num_taxa;
my @matrix;

parse_map_file();
process_whole_file();
print_file1s();
print_file2s();

#########################

sub parse_map_file {

#declare local array
my @map;
#declare local variables
my $i=0;

	open (IN1, "<", $ARGV[0]) || die ("Error: $!\n");
	@map = <IN1>;
	close IN1;
	while ($map[$i]){
		$line = $map[$i];
		chomp $line;
		@line = split (/\t/,$line);
		push (@species,$line[0]);
		push (@taxa_ids,$line[1]);
		$i++;
	}
}
print "\n@species\n"; #test
print "@taxa_ids\n\n\n\n\n";#test
#########################

sub process_whole_file {

#declare local variables
my $i=0;

	open (IN2, "<", $ARGV[1]) || die ("Error: $!\n");
	@whole = <IN2>;
	close IN2;
	while ($whole[$i]){
		$line = $whole[$i];
		chomp $line;
		if ($i==0){
			@line = split(/\s{1}/,$line);
			$num_taxa = $line[0];
			$num_char = $line[1];
		}
		elsif ($i==1){
			$num_species = $line;
		}
		elsif ($i==2){
			@species_ids = split(/\s{1}/,$line);
		}
		elsif ($i==3){
			@species_num_taxa = split(/\s{1}/,$line);
		}
		elsif ($i>3){
			push (@matrix,$line);
		}
		$i++;
	}
}
print "num_taxa $num_taxa\nnum_char $num_char\nnum_species $num_species\n";#test
print "@species_ids\n";#test
print "@species_num_taxa\n";#test
print "@matrix\n";#test
##########################

sub print_file1s {

#declare local variables
my $i=0;#index value
my $k;#adjust so file number matches taxon number
my $file1_name;
my $new_num_taxa;
my $element;
my $original_species_num_taxa;
my $new_species_num_taxa;
my $j=0;
my $taxa_id_to_remove;

	while ($i< $num_taxa){
		$k = $i+1;
		$file1_name = "file1_".$k.".txt";
		open (FILE1, ">>", $file1_name) || die ("Error: $!\n");
		$new_num_taxa = $num_taxa-1;
		print FILE1 "$new_num_taxa $num_char\n";
		print FILE1 "$num_species\n";
		print FILE1 "@species_ids\n";
		$element = $species[$i];#from map file
		$original_species_num_taxa = $species_num_taxa[$element];
		$new_species_num_taxa = $original_species_num_taxa-1;#from whole file
		$species_num_taxa[$element]=$new_species_num_taxa;
		print FILE1 "@species_num_taxa\n";

		while ($matrix[$j]){
			$line = $matrix[$j];
			$taxa_id_to_remove = $taxa_ids[$i];
			unless ($line =~ /$taxa_id_to_remove/){
				print FILE1 "$line\n";
			}
			$j++;
		}
		close FILE1;
		$i++;
		$j=0;
		$species_num_taxa[$element] = $original_species_num_taxa;
	}
}

#############################

sub print_file2s {

#declare local variables
my $i=0;
my $j=0;
my $k;
my $file2_name;
my $taxa_id_to_use;

	while ($i < $num_taxa){
		$k = $i+1;
		$file2_name = "file2_".$k.".txt";
		open (FILE2,">>",$file2_name) || die ("Error: $!\n");

		while ($matrix[$j]){
			$line = $matrix[$j];
			$taxa_id_to_use = $taxa_ids[$i];
			if ($line =~ /$taxa_id_to_use/){
				print FILE2 "$line\n";
			}
			$j++;
		}
		close FILE2;
		$i++;
		$j=0;
	}
}
