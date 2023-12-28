#!/usr/bin/perl
#Terri Porter, Aug.27, 2010
#Script to compile queryID, hitID, and cluster size from best_hit_summary (local database search, ublast) and results_parsed (dereplication)
#usage $perl hit_read_freq.plx top_blast.hits cluster.table
#modified to search for accessions beginning with "F" instead of "G"
#
use strict;
use warnings;

#declare variables
my $line;
my $i=0; #indexing
my $current_id;
my $id;
my $ref;
my $cluster_size;
my $current_ref;
my $k; #tempvar
my $j=0;#indexing
my $total_frequency;

#declare arrays
my @line;
my @ids;
my @refs;
my @refs_nonredundant;
my @cluster_sizes;

#declare hash
my %refs;

get_list_of_ids_and_refs();
#dereplicate_array_of_refs();
get_list_of_refs_and_cluster_sizes();
dereplicate_array_of_refs();
print_hit_read_frequency_table();
unlink("temp.txt");

#subroutine to get filtered list of ids and refs

sub get_list_of_ids_and_refs {

open (IN1, '<', $ARGV[0]) || die ("Error: $!\n");

while (<IN1>){
	$line = $_;
	chomp $line;
	if ($line =~ /^(F|G)/) {
		@line = split (/\t/, $line);
		$id = $line[0];
		push (@ids, $id);
		$ref = $line[1];
		push (@refs, $ref);
	}
	else {
		next;
	}
}
close IN1;

}

#subroutine to dereplicate @refs

sub dereplicate_array_of_refs {

%refs = map {$_,1} @refs;
@refs_nonredundant = keys %refs;

}

#subroutine to get list of refs and cluster_sizes

sub get_list_of_refs_and_cluster_sizes {

open (TMP, '>>', "temp.txt") || die ("Error: $!\n");

while ($ids[$i]){
	$current_id = $ids[$i];
	open (IN2, '<', $ARGV[1]) || die ("Error: $!\n");
	while (<IN2>) {
		$line = $_;
		chomp $line;
		if ($line =~ /$current_id/) {
			@line = split (/\t/, $line);
			$cluster_size = $line[1];
			print TMP "$refs[$i]\t$cluster_size\n";
		}
		else {
			next;
		}
	}
	close IN2;
	$i++;
}
close TMP;

}

#subroutine to print hit read frequency table

sub print_hit_read_frequency_table {

open (OUT, '>>', "hit_read_freq.txt") || die ("Error: $!\n");
print OUT "ReferenceID\tReadFrequency\n";

while ($refs_nonredundant[$j]){
	$current_ref = $refs_nonredundant[$j];
	open (IN3, '<', "temp.txt") || die ("Error: $!\n");
	while (<IN3>) {
		$line = $_;
		chomp $line;
		if ($line =~ /$current_ref/) {
			@line = split (/\t/,$line);
			$cluster_size = $line[1];
			push (@cluster_sizes,$cluster_size);
		}
		else {
			next;
		}
	}
	$j++;
	$total_frequency =0;
	foreach $k (@cluster_sizes) {
		$total_frequency += $k;
	}
	print OUT "$current_ref\t$total_frequency\n";
	@cluster_sizes=();
}
close IN3;
close OUT;

}
