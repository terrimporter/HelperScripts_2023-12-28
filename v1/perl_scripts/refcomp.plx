#!/usr/bin/perl
# Teresita M. Porter, June 8, 2019

# Script to compare three testNBC.fasta files by taxon (phylum & Arthropoda class) to create a table for R
# Create heatmap of taxa x refset filled with num seq in R
# USAGE perl refcomp.plx testNBC.fasta.NCBI testNBC.fasta.BOLD testNBC.fasta.v4

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $count;
my $phylum;
my $class;

# declar array
my @ncbi;
my @bold;
my @v4;
my @line;

# declare hash
my %ncbi_phyla; # key = phyla, value = seq count
my %ncbi_arthclass; # key = arthclass, value = seq count
my %bold_phyla; # key = phyla, value = seq count
my %bold_arthclass; # key = arthclass, value = seq count
my %v4_phyla; # key = phyla, value = seq count
my %v4_arthclass; # key = arthclass = value = seq count

open (IN1, "<", $ARGV[0]) || die "Cannot open first infile: $!\n";
@ncbi = <IN1>;
close IN1;

# hash infile
while ($ncbi[$i]) {
	$line = $ncbi[$i];
	chomp $line;

	if ($line =~ /^>/) { #header
		@line = split(/;/,$line);
		$phylum = $line[3];
		$class = $line[4];
		if (exists $ncbi_phyla{$phylum}) {
			$count = $ncbi_phyla{$phylum};
			$ncbi_phyla{$phylum} = $count + 1;
		}
		else {
			$ncbi_phyla{$phylum} = 1;
		}

		if ($phylum eq "Arthropoda") {
			if (exists $ncbi_arthclass{$class}) {
				$count = $ncbi_arthclass{$class};
				$ncbi_arthclass{$class} = $count + 1;
			}
			else {
				$ncbi_arthclass{$class} = 1;
			}
		}
	}
	else {
		$i++;
		next;
	}
	$i+=2;
}
$i=0;

open (IN2, "<", $ARGV[1]) || die "Cannot open second infile: $!\n";
@bold = <IN2>;
close IN2;

# hash infile
while ($bold[$i]) {
	$line = $bold[$i];
	chomp $line;

	if ($line =~ /^>/) { #header
		@line = split(/;/,$line);
		$phylum = $line[3];
		$class = $line[4];
		if (exists $bold_phyla{$phylum}) {
			$count = $bold_phyla{$phylum};
			$bold_phyla{$phylum} = $count + 1;
		}
		else {
			$bold_phyla{$phylum} = 1;
		}

		if ($phylum eq "Arthropoda") {
			if (exists $bold_arthclass{$class}) {
				$count = $bold_arthclass{$class};
				$bold_arthclass{$class} = $count + 1;
			}
			else {
				$bold_arthclass{$class} = 1;
			}
		}
	}
	else {
		$i++;
		next;
	}
	$i+=2;
}
$i=0;

open (IN3, "<", $ARGV[2]) || die "Cannot open third infile: $!\n";
@v4 = <IN3>;
close IN3;

# hash infile
while ($v4[$i]) {
	$line = $v4[$i];
	chomp $line;

	if ($line =~ /^>/) { #header
		@line = split(/;/,$line);
		$phylum = $line[3];
		$class = $line[4];
		if (exists $v4_phyla{$phylum}) {
			$count = $v4_phyla{$phylum};
			$v4_phyla{$phylum} = $count + 1;
		}
		else {
			$v4_phyla{$phylum} = 1;
		}

		if ($phylum eq "Arthropoda") {
			if (exists $v4_arthclass{$class}) {
				$count = $v4_arthclass{$class};
				$v4_arthclass{$class} = $count + 1;
			}
			else {
				$v4_arthclass{$class} = 1;
			}
		}
	}
	else {
		$i++;
		next;
	}
	$i+=2;
}
$i=0;

# print out a table in the right format for ggplot [R]

open (OUT, ">>", "refcomp.csv") || die "Error cannot open outfile :$!\n";
# refset, rank, taxon, count

# parse through NCBI phylum hash
while ( ($phylum, $count) = each (%ncbi_phyla)) {
	print OUT "NCBI,phylum,".$phylum.",".$count."\n";
}

# parse through NCBI arthclass hash
while ( ($class, $count) = each (%ncbi_arthclass)) {
	print OUT "NCBI,arthclass,".$class.",".$count."\n";
}

# parse through BOLD phylum hash
while ( ($phylum, $count) = each (%bold_phyla)) {
	print OUT "BOLD,phylum,".$phylum.",".$count."\n";
}

# parse through BOLD arthclass hash
while ( ($class, $count) = each (%bold_arthclass)) {
	print OUT "BOLD,arthclass,".$class.",".$count."\n";
}

# parse through v4 phylum hash
while ( ($phylum, $count) = each (%v4_phyla)) {
	print OUT "v4,phylum,".$phylum.",".$count."\n";
}

# parse through v4 arthclass hash
while ( ($class, $count) = each (%v4_arthclass)) {
	print OUT "v4,arthclass,".$class.",".$count."\n";
}

close OUT;
