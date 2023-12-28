#!/usr/bin/perl
# Teresita M. Porter, Feb. 3/20
# Script to parse fasta header, grab taxon name, check for any non-nucleotide ambiguities
# USAGE perl parse_pseudogenes.plx pseudogenes.fasta.strict bold.nt.fasta

use strict;
use warnings;
use Data::Dumper;

# var
my $i=0;
my $line;
my $accession;
my $genus;
my $species;
my $binomialname;
my $j;
my $seq;
my $lineage;
my $oldrecord;
my $outfile;
my $records;
my $record;

# array
my @in;
my @line;
my @in2;
my @lineage;
my @records;
my @split;

# hash
my %pseudogene; # key = binomial name, value = accession\n seq
my %bold; # key = species, value = accession lineage\n seq

open (IN, "<", $ARGV[0]) || die "Error cannot open infile:$!\n";
@in = <IN>;
close IN;

# hash pseudogene name and seq
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;

		unless ($line =~ /aff\.|sp\.|cf\./) { # only work with fully identified pseudogenes
			@line = split(/ /, $line);
			$accession = $line[0];
			$genus = $line[1];
			$species = $line[2];
			$binomialname = $genus." ".$species;
			$j = $i+1;
			$seq = $in[$j];
			chomp $seq;

			if (exists $pseudogene{$binomialname}) {
				$oldrecord = $pseudogene{$binomialname};
				$pseudogene{$binomialname} = $oldrecord."|".$accession."\n".$seq;
			}
			else {
				$pseudogene{$binomialname} = $accession."\n".$seq;
			}
		}

	}
	$i++;
}
$i=0;

# print Dumper(\%pseudogene);

open (IN2, "<", $ARGV[1]) || die "Error cannot open infile2: $!\n";
@in2 = <IN2>;
close IN2;

# parse the bold file and hash the lineage, name, and seq
while ($in2[$i]) {
	$line = $in2[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;

		unless ($line =~ /aff\.|sp\.|cf\./) {
			@line = split(/ /,$line);
			$accession = shift @line; # remove accession
			$lineage = join ' ', @line;
			@lineage = split(/;/, $lineage);
			$species = pop @lineage;

			if (exists $pseudogene{$species}) { # if find pseudogene species in BOLD
				$j = $i+1;
				$seq = $in2[$j];
				chomp $seq;
				if (exists $bold{$species}) {
					$oldrecord = $bold{$species};
					$bold{$species} = $oldrecord."|".$accession." ".$lineage.";".$species."\n".$seq;
				}
				else {
					$bold{$species} = $accession." ".$lineage.";".$species."\n".$seq;
				}
			}
		}
	}
	$i++;
}
$i=0;

# print Dumper(\%bold);

# for each species in %bold create new file and print pseudogene and BOLD seqs
while ( ($species, $record) = each %bold) {
	$species =~ s/ /_/g;
	$outfile = $species.".fasta";
	open (OUT, ">>", $outfile) || die "Cannot open outfile for species $species\n";

	$species =~ s/_/ /g;
	$records = $bold{$species};
	@records = split(/\|/, $records);

	foreach $record (@records) {
		print OUT ">".$record."\n";
	}

	$records = $pseudogene{$species};
	@records = split(/\|/, $records);

	foreach $record (@records) {
		@split = split(/\n/,$record);
		$accession = $split[0];
		$seq = $split[1];

		print OUT ">".$accession." pseudogene\n".$seq."\n";
	}
	close OUT;
}
