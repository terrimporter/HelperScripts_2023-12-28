#!/usr/bin/perl
#Nov.16,2010 by Terri Porter
#Script to sort file.qual by the IDs in the sorted .fna and remove first 10 bases
#usage $perl sort_qual_by_agrepID.plx file.qual sorted_trimmed.fna
#modified Nov.18,2010 to sort qual by primer sorted fasta, but not to trim.  Use with MIDX_trimmed.qual
#modified Nov.24,2010 to automatically remove new line endings from within the sequence part
#modify again to do order_fasta_to_qual.plx stuff too

use strict;
use warnings;

#declare var
my $line;
my $id_to_match;
my $id;
my $x;
#my $match=0;
my $i=1;
my $trimmed_phred;
my $flag=0;
my $w=0;
my $phred;
my $header;

#declare array
my @qual;
my @fna;
my @ids;
my @line;

open (IN1,"<",$ARGV[0]) || die ("Error: $!\n");
open (TEMP1,">>","temp1.txt") || die ("Error: $!\n");

#remove new line character from sequence-part
while (<IN1>) {
	$line =$_;
	chomp $line;
	if ($line =~ /^>/) {
		print TEMP1 "\n",$line,"\n";
		$flag=0;
	}
	elsif ($flag==0) {
		print TEMP1 $line;
		$flag=1;
	}
	elsif ($flag>0) {
		print TEMP1 " $line";
	}
}
close IN1;
close TEMP1;

#remove empty line at top of file from before
open (IN2,"<","temp1.txt") || die ("Error: $!\n");
open (TEMP2,">>","temp2.txt") || die ("Error: $!\n");

while (<IN2>) {
	$line = $_;
	chomp $line;
	if ($line =~ /\S+/) {
		print TEMP2 $line,"\n";
	}
}
close IN2;
close TEMP2;

#grab qual entries based on sorted_dereplciated.fasta
open (QUAL,"<","temp2.txt") || die ("Error: $!\n");
@qual = <QUAL>;

open (FNA,"<",$ARGV[1]) || die ("Error: $!\n");
while (<FNA>) {
	$line = $_;
	chomp $line;
	if ($line =~ />/) {
		$line =~ />(\w{14})/;
		$id_to_match = $1;
		push(@ids, $id_to_match);
	}
}
close FNA;

open (OUT,">>","sorted.qual") || die ("Error: $!\n");

while ($qual[$w]) {
	$line = $qual[$w];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w{14})/;
		$id = $1;
		foreach $id_to_match (@ids) {
			if ($id_to_match eq $id) {
				$x = $w+1;
				$phred = $qual[$x];
				print OUT ">$id\n$phred\n";
			}
		}
	}
	$w++;
}
close QUAL;
close OUT;
unlink("temp1.txt","temp2.txt");

order_fasta_to_qual();

sub order_fasta_to_qual {

	#var
	my $i=0;
	my $line;
	my $id;
	my $j=0;
	my $line2;
	my $id2;
	my $k;
	my $seq;

	#array
	my @fasta;
	my @qual;

	open (FASTA,"<", $ARGV[1]) || die ("Error: $!\n");
	@fasta = <FASTA>;
	close FASTA;

	open (QUAL2,"<","sorted.qual") || die ("Error: $!\n");
	@qual=<QUAL2>;
	close QUAL2;

	open (OUT2,">>","sorted.fasta") || die ("Error: $!\n");

	while ($qual[$i]) {
		$line = $qual[$i];
		chomp $line;
		if ($line =~ /^>/) {
			$line =~ /^>(\w{14})\s*/;
			$id = $1;
			$j=0;
			my $flag=0;
			while ($fasta[$j]) {
				$line2 = $fasta[$j];
				chomp $line2;
				if ($line2 =~ /^>/) {
					$line2 =~ /^>(\w{14})\s*/;
					$id2 = $1;
					if ($id2 eq $id && $flag==0) {
						$k = $j+1;
						$seq = $fasta[$k];
						chomp $seq;
						print OUT2 ">$id2\n$seq\n";
						$flag=1;
					}
					elsif ($id2 eq $id && $flag==1) {
						$j++;
						next;
					}
				}
				$j++;
			}
		}
		$i++;
	}
}
