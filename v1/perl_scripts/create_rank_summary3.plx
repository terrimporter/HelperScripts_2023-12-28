#!/usr/bin/perl
#March 31, 2011 by Terri Porter
#Script to take gi\tgb\tlineage\n file and create gi_original\trank\n file
#usage $perl create_rank_summary.plx hit_lineage.txt gb_gi.query ITS.blast.parsed.all
#modified May 11, 2011 to work with hit_lineage.txt.mapped format query_gi|query_gb\thit_gb\tlineage
#usage $perl create_rank_summary2.plx query_lineage.txt hit_lineage.txt.mapped
#fix so that it works right!!!!bug - if no family, then no flag=1, then don't grab genus

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $j=0;
my $k=0;
my $lineage;
my $rank;
my $l=0;
my $genus_maybe;
my $genus="nil";
my $family="nil";
my $order="nil";
my $class="nil";
my $phylum="nil";
my $kingdom="nil";
my $hit_gi;
my $hit_gb;
my $index;
my $filename;
my $next=0;
my $gi_part;
my $z=0;
my $end=0;
my $m=0;
my $map_gi;
my $query_gi;
my $query_gb;
my $line2;
my $temp;
my $species="nil";
my $species_maybe;
my $query;
my $rank_name;
my $query_lineage;
my $counter=0;
my $mapping;
my $flag=0;
my $num;
my $index1;
my $undef;

#declare array
my @hit_lineage;
my @blast;
my @original_gb;
my @hit_full_gb;
my @line;
my @lineage;
my @lin;
my @gi_part;
my @map;
my @mapped_gis;
my @query;
my @hit;
my @temp;
my @line2;

open (QUERY,"<",$ARGV[0]) || die ("Error cannot read from query_lineage.txt: $!\n");
@query = <QUERY>;
close QUERY;

open (HIT,"<",$ARGV[1]) || die ("Error cannot read from hit_lineage.txt.mapped: $!\n");
@hit = <HIT>;
close HIT;

open (TEMP,">>","hit_lineage.parsed") || die ("Error cannot write to hit_lineage.temp: $!\n");

while ($query[$i]) {
	$line = $query[$i];
	chomp $line;
	$line =~ /^(\d+)/;
	$query_gi = $1;
	while ($hit[$j]) {
		$line2 = $hit[$j];
		chomp $line2;
		@line2 = split(";",$line2);
		$temp = shift(@line2);
		@temp = split(/\t/,$temp);
		pop(@temp);#remove Eukaryota;
		$mapping = join("\t",@temp);
		$num = scalar(@line2);
		$index1 = $num-1;
		$undef = $num+1;

		if ($line2 =~ /^$query_gi/){
			#parse through lineage starting with species
			while ($line2[$index1]) {
			if ($index1>=0) {
				$rank = $line2[$index1];
				#print "$rank\n";#test
				$rank =~ s/^\s+//; 
				$rank =~ s/incertae sedis//;
				
				if ($next==1) {
					$rank =~ s/\;//;
					$genus_maybe = $rank;
						
					if ($genus_maybe =~ /^[A-Z]/) {
						$genus = $genus_maybe;
						#		print $genus."\n";
					}
					$next=0;
				}
				elsif ($next==0) {
					if ($rank =~ /\./) {
						$species_maybe = $rank;

						if ($species_maybe =~ /^[a-z]/) {
							$species = $species_maybe;
							$species =~ s/\.//;
							#	print $species."\n";
						}
						$next=1;
					}
					elsif ($rank =~ /aceae$/) {
						$family = $rank;
						#	print $family."\n";
					}
					elsif ($rank =~ /ales$/) {
						$order = $rank;
					}
					elsif ($rank =~ /Ichthyophonida/) {
						$order = $rank;
					}
					elsif ($rank =~ /mycetes$/) {
						$class = $rank;
					}
					elsif ($rank =~ /Ichthyosporea/) {
						$class = $rank;
					}
					elsif ($rank =~ /mycota$/) {
						$phylum = $rank;
					}
					elsif ($rank =~ /Fungi$/) {
						$kingdom = $rank;
					}
				}
				$index1--;
			}
			#print "completed if statement\n";
			else {
				$index1=$undef;
			}
			}
			#print "completed while loop\n";
			#$z=0;
			#print "almost to flag zero\n";#test
			if ($flag==0) {
				#	print "got to flag zero\n";#test
				print TEMP "$mapping\t$species\t$genus\t$family\t$order\t$class\t$phylum\t$kingdom\n";
				$species="nil";$genus="nil";$family="nil";$order="nil";$class="nil";$phylum="nil";$kingdom="nil";
				$next=0;
				$flag=1;
			}
		}	
		$j++;
	}
	$j=0;
	@line2=();@temp=();$line2=();$temp=();$mapping=();
	$i++;
	$flag=0;
}
close TEMP;

print "Please enter what rank to summarize (ex. species, genus, family, order, class, phylum, kingdom):\n";
$rank = <STDIN>;
chomp $rank;

if ($rank eq 'species') {
	$index = 3;
}
if ($rank eq 'genus') {
	$index = 4;
}
elsif ($rank eq 'family') {
	$index=5;
}
elsif ($rank eq 'order') {
	$index=6;
}
elsif ($rank eq 'class') {
	$index=7;
}
elsif ($rank eq 'phylum') {
	$index=8;
}
elsif ($rank eq 'kingdom') {
	$index=9;
}

open (LIN,"<","hit_lineage.parsed") || die ("Error cannot read from hit_lineage.temp: $!\n");
@lin = <LIN>;
close LIN;

$j=0;

while ($query[$j]) {
	$query = $query[$j];
	chomp $query;
	if ($query =~ /^\d+/) {
		$query =~ /^(\d+)/;
		$query_gi = $1;
#		print $query_gi."\n";
	}

	while ($lin[$k]) {
		$line = $lin[$k];
		chomp $line;
#		print $line."\n";		
		if ($line =~ /^$query_gi/) {
#			print "match";#test
			@line = split(/\t/,$line);
			$rank_name = $line[$index];
			print $rank_name."\n";#test
			if ($query =~ /$rank_name/){
				$counter++;
			}
		}
		$k++;
	}
	$k=0;
	$line=();
	$query=();
	$j++;
}
print "$rank rank correct: $counter\n";
