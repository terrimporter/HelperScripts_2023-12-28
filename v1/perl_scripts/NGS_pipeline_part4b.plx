#!/usr/bin/perl
#Feb.8,2011 by Terri Porter
#Script to parse -outfmt 0 (pairwise text) and grab hit_gb for fully identified top hit
#Modified May 12, 2011 usage $perl NGS_pipeline_part4c.plx

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

#var
my $first;
my $outfile2b;
my $i=0;
my $line;
my $flag=0;
my $query_map;
my $hit_gb;
my $j=0;
my $factory;
my $k=0;
my $hit_gb;
my $hit_gi;
my $species;
my $lineage;
my $l=0;
my $m=0;
my $gb;
my $flag2=0;

#array
my @blast;
my @query_map;
my @line;
my @in2;
my @gb;
my @hit_gis;
my @read;
my @organism;
my @read2;

#parse blast+ -outfmt 0 blast results to grab hit_gb for a fully identified top hit
print "Enter name of file.blast:\n";
$first = <STDIN>;
chomp $first;

open (IN,"<",$first) || die ("Error cannot read blast file: $!\n");
@blast = <IN>;
close IN;

$outfile2b = "query_map_hit_gb.txt";

open (OUTb,">>",$outfile2b) || die ("Error: $!\n");

while ($blast[$i]) {
	$line = $blast[$i];
	chomp $line;
	if ($flag==0) {
		if ($line =~ /^Query\=/) {
			$line =~ /Query\=\s+(\d+\|\w+)/g;###edit regex here###
			$query_map = $1;
			push(@query_map,$query_map);
			$flag=1;
			print "Processing query: $query_map\n";#test
			print OUTb "$query_map\t";
		}	
	}
	elsif ($flag==1) {
		if ($line =~ /^\w+\|\w+\.\d+\|/g) {
			if ($line =~ /(Uncultured|uncultured|endophyte|mycorrhizae|mycorrhizal samples|unclassified|symbionts)/) {
				$i++;
				next;
			}
			elsif ($line =~ /sp\./) {
				#print "match\n";#test
				$i++;
				next;
			}
			elsif ($line =~ /cf\./) {
				$i++;
				next;
			}
			elsif ($line =~ /aff\./) {
				$i++;
				next;
			}
			elsif ($line =~ /environmental sample/) {
				$i++;
				next;
			}
			else {
				@line = split(/\|/,$line);
				$hit_gb = $line[1];
				print OUTb "$hit_gb\n";
				$flag=0;
			}
			@line=();
		}
	}
	$i++;
}
close OUTb;

#grab hit_gb and grab genbank files for each (.txt) then parse out lineage

open (IN2,"<","query_map_hit_gb.txt") || die ("Error cannot read query_map_hit_gb.txt: $!\n");
@in2 = <IN2>;
close IN2;

while ($in2[$j]) {
	$line = $in2[$j];
	chomp $line;
	@line = split(/\t/,$line);
	$query_map = $line[0];
	push(@query_map,$query_map);
	$hit_gb = $line[1];
	push (@gb,$hit_gb);
	$j++;
}

$factory = Bio::DB::EUtilities-> new(	-eutil	=> 'efetch',
					-db	=>'nucleotide',
					-id	=> \@gb,
					-rettype	=> 'gb');

@hit_gis = split(m{\n},$factory->get_Response->content);

open (GB_TEXT,">>","hit_gb.txt") || die ("Error cannot write to hit_gb.txt: $!\n");

foreach $line (@hit_gis) {
	print GB_TEXT $line."\n";
}
close GB_TEXT;

open (IN2,"<","hit_gb.txt") || die ("Error cannot read from hit_gb.txt: $!\n");
@read = <IN2>;
close IN2;

open (OUT2,">>","hit_lineage.txt") || die ("Error cannot write to hit_lineage.txt:$!\n");

$flag=0;

while ($read[$k]) {
	$line = $read[$k];
	chomp $line;

	if ($flag==0) {
		if ($line =~ /VERSION/) {
			$line =~ /VERSION\s+(\w+\.\d+)\s+GI:(\d+)/;
			$hit_gb = $1;
			$hit_gi = $2;
		}
		elsif ($line =~ /\s{2}ORGANISM/) {
			$line =~ s/  ORGANISM  //;
			@line = split(/ /,$line);
			$species = $line[1];
			$flag=1;
		}
	}
	elsif ($flag==1) {
		if ($line !~ /(REFERENCE|COMMENT)/) {
			push(@organism, $line);
		}
		elsif ($line =~ /(REFERENCE|COMMENT)/) {
			$lineage = join(" ",@organism);
			$lineage =~ s/            //g;
			$lineage =~ s/.$/;/;
			$flag=0;
			print OUT2 "$hit_gi\t$hit_gb\t$lineage. $species..\n";#test
			$hit_gb=();
			$hit_gi=();
			@organism=();
			$lineage=();
		}
	}
	$k++;
}
close OUT2;

open (IN3,"<","hit_lineage.txt") || die ("Error cannot read from hit_lienage.txt: $!\n");
@read2 = <IN3>;
close IN3;

open (OUT3,">>","hit_lineage.txt.mapped") || die ("Error cannot write to hit_lineage.txt.mapped: $!\n");

while ($read2[$l]) {
	$line = $read2[$l];
	chomp $line;

	@line = split(/\t/,$line);
	$hit_gb = $line[1];

	while ($gb[$m]) {
		$gb = $gb[$m];
		
		if ($gb =~ /$hit_gb$/) {
			$query_map = $query_map[$m];
			
			if ($flag2 == 0) {
				print OUT3 "$query_map\t$line\n";	
				$flag2=1;
				$gb[$m]="nil";
			}
		}
		$m++;
	}
	$m=0;
	$l++;
	$flag2=0;
}
