#!/usr/bin/perl
#Terri Porter, Jan. 15, 2013
#Script to parse out Spizellomyces punctatus sequences from ortholog alignments from trimAl
#usage perl parse_out_Sp.plx file.fasta
#be sure that file.fasta contains sequence on a single line

use strict;
use warnings;

#declare var
my $line;
my $i=0;
my $infile;
my $orthoId;
my $seqId;
my $flag=0;
my $seq='';

#declare array
my @in;
my @infile;
my @orthoIds;

#declare hash
my %SP_orthoId_seqId;
my %AB_orthoId_seqId;
my %SP_orthoId_seq;
my %AB_orthoId_seq;

open (IN, "<", $ARGV[0]) || die "Error cannot read infile: $!\n";
@in = <IN>;
close IN;

#parse orthoid from filename
$infile = $ARGV[0];
chomp $infile;
if ($infile !~ /^\./ && $infile =~ /^\d+/) {
	@infile = split(/\./, $infile);
	$orthoId = $infile[0];
	push(@orthoIds, $orthoId);
}

#parse seqid and seq
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>(SP|AB)/) { # SP - 454, AB - Broad
		$line =~ /^>(\S+)/;
		$seqId = $1;
		if ($seqId =~ /^SP/) {
			$SP_orthoId_seqId{$orthoId} = $seqId;
			#print "seqid $seqId\n";
		}
		elsif ($seqId =~ /^AB/) {
			$AB_orthoId_seqId{$orthoId} = $seqId;
		}

		#$flag=1;
	}
	else {
		#if ($seq eq '') { # not defined
		$seq = $line;
		if ($seqId =~ /^SP/) {
			$SP_orthoId_seq{$orthoId} = $seq;
		}
		elsif ($seqId =~ /^AB/) {
			$AB_orthoId_seq{$orthoId} = $seq;
		}

		#}
		#elsif ($seq ne '') { #defined
		#$seq = $seq.$line; #account for sequence on multiple lines
		#}
	}
#	elsif ($flag==1 && $line =~ /^>/) { #hash
#		$flag=0;
#		if ($seqId =~ /^SP/) {
#			$SP_orthoId_seq{$orthoId} = $seq;
#		}
#		elsif ($seqId =~ /^AB/) {
#			$AB_orthoId_seq{$orthoId} = $seq;
#		}
		#$orthoId=();
#		$seq='';
	
#		if ($line =~ /^>(SP|AB)/) {
#			$line =~ /^>(\S+)/;
#			$seqId = $1;
#			if ($seqId =~ /^SP/) {
#				$SP_orthoId_seqId{$orthoId} = $seqId;
#			}
#			elsif ($seqId =~ /^AB/) {
#				$AB_orthoId_seqId{$orthoId} = $seqId;
#			}	
#			$flag=1;
#		}
#	}
	$i++;
}
$i=0;

#print presence absence table
open (OUT1, ">>", "presence_absence.txt") || die "Error cannot open outfile1: $!\n";

while($orthoIds[$i]) {
	$orthoId = $orthoIds[$i];
	chomp $orthoId;

	print OUT1 $orthoId."\t";

	if (exists $SP_orthoId_seqId{$orthoId}) {
		print OUT1 "1\t";
	}
	else {
		print OUT1 "0\t";
	}

	if (exists $AB_orthoId_seqId{$orthoId}) {
		print OUT1 "1\n";
	}
	else {
		print OUT1 "0\n";
	}
	$i++;
}
$i=0;
close OUT1;

#print orthoid seqid map
#open (OUT2, ">>", "orthoid_seqid.map") || die "Error cannot open outfile2: $!\n";

#while($orthoIds[$i]) {
#	$orthoId = $orthoIds[$i];
#	chomp $orthoId;

#	print OUT2 $orthoId."\t";

#	if (exists $SP_orthoId_seqId{$orthoId}) {
#		$seqId = $SP_orthoId_seqId{$orthoId};
#		print OUT2 "$seqId\t";
#	}
#	else {
#		print OUT2 "\t";
#	}

#	if (exists $AB_orthoId_seqId{$orthoId}) {
#		$seqId = $AB_orthoId_seqId{$orthoId};
#		print OUT2 "$seqId\n";
#	}
#	else {
#		print OUT2 "\n";
#	}
	
#	$i++;
#}
#$i=0;
#close OUT2;

#print fasta file for blasting >$seqId\n$seq\n
open (OUT3, ">>", "SP.fasta") || die "Error cannot open outfile3: $!\n";

while (($orthoId,$seq) = each (%SP_orthoId_seq)) {
	if (exists $SP_orthoId_seq{$orthoId}) {
		print OUT3 ">$orthoId\n$seq\n";
	}
	else {
		next;
	}
}
close OUT3;

open (OUT4, ">>", "AB.fasta") || die "Error cannot open outfile4: $!\n";

while(($orthoId,$seq) = each (%AB_orthoId_seq)) {
	if (exists $AB_orthoId_seq{$orthoId}) {
		print OUT4 ">$orthoId\n$seq\n";
	}
	else {
		next;
	}
}
close OUT4;






