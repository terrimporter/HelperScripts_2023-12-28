#!/usr/bin/perl
#August 23, 2012 by Terri Porter
#Script to parse trimal.nex to figure out % char (num char / length ortho) for each taxon to create a heatmap
#usage perl parse_nexus.plx trimal.nex

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $flag=0;
my $ortho;
my $taxon;
my $seq;
my $orthoLength;
my $char;
my $j=0;
my $count=0;
my $prop; #proportion=(char/orthoLength)*100

#declare array
my @in;
my @line;
my @seq;

#declare hash
my %ortho_length;
my %ortho_taxon_count; #hash of hashes

open (IN, "<", $ARGV[0]) || die "Error cannot open trimal.nex: $!\n";
@in = <IN>;
close IN;

#parse number of char per ortho per taxon
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0) {
		if ($line =~ /\[\d+\]/) {
			$line =~ /\[(\d+)\]/;
			$ortho = $1;
			#print "$ortho\t";#test
			$flag=1;
		}
	}

	elsif ($flag==1) {
		if ($line !~ /^\S+/) {
			$flag=2;
			#print "got to flag 2\n";#test
		}
		#print "failed to recognize empty line\n";#test
	}

	elsif ($flag==2) {
		if ($line =~ /^[A-Z][A-Z]\s+/) {
			@line = split(/\t/,$line);
			$taxon = $line[0];
			#print "$taxon\t";#test
			$seq = $line[1];
			@seq = split (//,$seq);
			$orthoLength = scalar(@seq);
			#print "$orthoLength\t";#test
			$ortho_length{$ortho} = $orthoLength;

			while($seq[$j]) {
				$char = $seq[$j];
				if ($char =~ /\-/) {
					$j++;
					next;
				}
				else {
					$count++;
				}
				$j++;
			}
			$j=0;
			$char=();
			$ortho_taxon_count{$ortho}{$taxon} = $count;
			#print "$count\n";#test
			$count=0;
		}

		elsif ($line !~ /^\S+/) {
			$flag=0;
			#	print "reset flag to 0\n";#test
		}
		#print "did not recognize empty line after alignment\n";#test
	}

	$i++;

	$line=();
	@line=();
	$taxon=();
	$seq=();
	@seq=();
	$orthoLength=();
}
$i=0;

#print file for R

open (OUT, ">>", "heatmap.txt") || die "Error cannot open outfile: $!\n";

print OUT "\tAA\tAB\tAM\tAN\tBB\tBD\tBE\tCA\tCR\tEP\tGP\tLB\tMB\tMC\tMV\tNE\tPB\tPD\tPE\tRO\tSC\tSP\tSR\tUM\n";

while ( ($ortho, $taxon) = each %ortho_taxon_count) {

	print OUT "$ortho\t";

	$count = $ortho_taxon_count{$ortho}{'AA'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'AB'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'AM'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'AN'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'BB'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'BD'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'BE'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'CA'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'CR'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'EP'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'GP'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'LB'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'MB'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'MC'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'MV'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'NE'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'PB'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'PD'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'PE'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'RO'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'SC'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'SP'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'SR'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\t";

	$count = $ortho_taxon_count{$ortho}{'UM'};
	$orthoLength = $ortho_length{$ortho};
	$prop = ($count/$orthoLength)*100;
	printf OUT ("%.1f", $prop);
	print OUT "\n";

}

close OUT;
