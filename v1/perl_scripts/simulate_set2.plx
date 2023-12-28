#!/usr/bin/perl
#Sept.23,2011 by Terri Porter
#Script to take features.txt.parsed from parse_features_with_gb_list.pxl and a file with primer sequences (format: primername\tprimersequence(5'-3')\tFWD/REV\n)
#find primer sites, grab fragments of various sizes (50, 100, 200, and 400)
#usage perl simulate_set.plx primers.txt features.txt.parsed

use strict;
use warnings;
use Text::Levenshtein qw(fastdistance);

#declare var
my $primerlength;

#declare array
my @primer;
my @features;
my @primersequence;
my @primer_mismatch;

#declare hash
my %primer;
my %orientation;
my %organism;
my %classification;
my %LSU_seq;

open (PRIMER, "<", $ARGV[0]) || die "Cannot read primers.txt: $!\n";
@primer = <PRIMER>;
close PRIMER;

process_primers();

open (FEATURES,"<",$ARGV[1]) || die "Cannot read features.txt.parsed: $!\n";
@features = <FEATURES>;
close FEATURES;

process_features();

find_primers();

#remove_duplicates();

#####################

sub process_primers{
print "Processing primers...\n"; #status
#declare var
my $a=0;
my $line;
my $primername;
my $primersequence;
my $orientation;

#declare array
my @line;

while ($primer[$a]) {
	$line = $primer[$a];
	chomp $line;
	@line = split(/\t/,$line);
	$primername = $line[0];
	$primersequence = $line[1];
	$orientation = $line[2];
	$primer{$primername} = $primersequence;
	$orientation{$primername} = $orientation;
	$a++;
}

}

####################

sub process_features {
print "Processing features...\n";#status
#declare var
my $b=0;
my $line;
my $accession;
my $organism;
my $classification;
my $LSU_seq;

#declare array
my @line;

open (FASTAOUT,">>","full.fasta") || die "Error cannot write to full.fasta :$!\n";

while ($features[$b]) {
	$line = $features[$b];
	chomp $line;
	@line = split(/\t/,$line);
	$accession = $line[0];
	$organism = $line[1];
	$classification = $line[5];
	$LSU_seq = $line[6];
	#print "LSU_seq:$LSU_seq\n";#test
	$organism{$accession} = $organism;
	$classification{$accession} = $classification;
	$LSU_seq{$accession} = $LSU_seq;
	print FASTAOUT ">$accession\n$LSU_seq\n";
	$b++;
}
close FASTAOUT;

}

####################

sub find_primers{

print "Finding primers...\n";#status

#declare var
my $primername;
my $primersequence;
#my $primerlength;
my $orientation;
my $filename;
my $accession;
my $LSU_seq;
my $LSU_seq_length;
my $c_temp;
my $c_max;
my $c=0;
my $window;
my $distance;
my $d=1;
my $base;
my $sub400;
my $c_rev;
my $f=0;
my $primer_mismatch;

#declare array
#my @primersequence;
my @LSU_seq;
my @substring;
my @substring_rev;

while (($primername,$primersequence) = each (%primer)) {
	@primersequence = split(//,$primersequence);
	$primerlength = scalar(@primersequence);
	$orientation = $orientation{$primername};
	print "\nprimer: $primername\n"; #status mesg
	#create_primer_mismatch();

	$filename = $primername."_400.fasta";
	open (OUT,">>",$filename) || die "Error cannot write to $filename:$!\n";
	
	while (($accession,$LSU_seq) = each (%LSU_seq)) {
		#print "$accession\n";#test
		#print ".\n";#status
		@LSU_seq = split(//,$LSU_seq);
		$LSU_seq_length = scalar(@LSU_seq);
		$c_max = $LSU_seq_length-$primerlength;
				
		while ($c <= $c_max) {
			#print ".";#status
			$window = substr $LSU_seq,$c,$primerlength;
			$distance = fastdistance($primersequence,$window);
			#$c_temp=$c;
			if ($distance <= 1) {###Adjust distance cutoff here###
				print "got fuzzy match\n";#status
				if ($orientation eq "FWD") {
					print "got a fwd one\n";#status
					while ($LSU_seq[$c]) {
						if ($d<= 400) {
							$base = $LSU_seq[$c];
							push(@substring,$base);
						}
						$c++;
						$d++;
					}
					$d=1;
					$sub400 = join("",@substring);
					print OUT ">$accession\n$sub400\n";
					$c=$c_max;
						
				}
				elsif ($orientation eq "REV") {
					print "got a rev one\n";#status
					$c_rev = $c+$primerlength-1;
					while ($LSU_seq[$c_rev]) {
						if ($d<=400) {
							$base = $LSU_seq[$c_rev];
							push(@substring,$base);
						}
						$c_rev--;
						$d++;
					}
					$d=1;
					@substring_rev = reverse(@substring);
					$sub400 = join("",@substring_rev);
					print OUT ">$accession\n$sub400\n";
					$c=$c_max;
				}
				@substring=();
				$sub400=();
				@substring_rev=();
			}
			#$c_temp=$c;
			#$f=0;
			$c++;
		}
		$c=0;
		#print ".";#test
	}	
	close OUT;
	@primersequence=();
	$primerlength=();
	$orientation=();
	@LSU_seq=();
	}

}
####################

sub create_primer_mismatch {

print "Creating primer mismatches...\n";#status

#declare var
my $e=0;
my $base_original;
my $x;
my $newprimer;

#declare array
my @newprimer;

@primer_mismatch=();

while ($e < $primerlength) {
	$base_original = $primersequence[$e];
	
	$primersequence[$e] = "A";
	foreach $x (@primersequence) {
		push(@newprimer,$x);
	}
	$newprimer = join("",@newprimer);
	push(@primer_mismatch,$newprimer);
	@newprimer=();
	
	$primersequence[$e] = "C";
		foreach $x (@primersequence) {
		push(@newprimer,$x);
	}
	$newprimer = join("",@newprimer);
	push(@primer_mismatch,$newprimer);
	@newprimer=();

	
	$primersequence[$e] = "T";
	foreach $x (@primersequence) {
		push(@newprimer,$x);
	}
	$newprimer = join("",@newprimer);
	push(@primer_mismatch,$newprimer);
	@newprimer=();


	$primersequence[$e] = "G";
	foreach $x (@primersequence) {
		push(@newprimer,$x);
	}
	$newprimer = join("",@newprimer);
	push(@primer_mismatch,$newprimer);
	@newprimer=();


	$primersequence[$e]=$base_original;
	$e++;
}
#print "\n";#test
}
####################

sub remove_duplicates {
print "Removing duplicates...\n";#status
#declare var
my $g=0;
my $filename;
my $h=0;
my $i;
my $j=0;
my $line;
my $accession;
my $fragment;
my $newfilename;
my $sequence;

#declare array
my @output;
my @infile;
my @unique;

#declare hash
my %fragment;

@output = qx(ls | grep "_400.fasta");

while ($output[$g]) {
	$filename = $output[$g];
	chomp $filename;
	print "Dereplicating $filename...\n";#status
	open (FILE,"<",$filename) || die "Error cannot open $filename:$!\n";
	@infile = <FILE>;
	close FILE;

	while ($infile[$h]) {
		$line = $infile[$h];
		chomp $line;

		if ($line=~ /^>/) {
			$line =~/^>(\w+)/;
			$accession = $1;
			$i=$h+1;
			$fragment = $infile[$i];
			chomp $fragment;
			$fragment{$accession} = $fragment;
			$h++;
		}
		$h++;
	}
	$h=0;
	@unique = keys %fragment;

	$newfilename = $filename.".nonredundant";
	open (OUT,">>",$newfilename) || die "Error cannot write to $newfilename:$!\n";

	while ($unique[$j]) {
		$accession = $unique[$j];
		$sequence = $fragment{$accession};
		print OUT ">$accession\n$sequence\n";
		$j++;
	}
	$j=0;
	close OUT;
	$g++;
}


}
