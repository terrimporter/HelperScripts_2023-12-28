#!/usr/bin/perl
#Oct. 5, 2011 by Terri Porter
#Script to parse gblocks.nex into phylip formatted files for each ortholog for ProtTEst
#usage perl nex_to_phylip.plx gblocks.nex

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $ntax;
my $nchar;
my $ortho_id;
my $taxon;
my $seq;
my $filename;
my $flag=0;
my $j=1;

#declare array
my @in;
my @line;
my @value;

#declare hash
my %ortho_aln;

open (IN,"<",$ARGV[0]) || die "Error cannot read from gblocks.nex: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	
	if ($flag==0) {
		if ($line =~ /^\[\d{8}\]/) {
			$line =~ /^\[(\d{8})\]/;
			$ortho_id = $1;

			$filename = "$ortho_id".".phy";
			open (OUT,">>",$filename) || die "Error cannot write to $filename:$!\n";
			#print OUT "$ntax\t$nchar\n";
			$flag=1;
		}
	}
	elsif ($flag>0) {
		#print "got here\n";#test
		if ($flag==1) {
			#print "got here\n";#test
			if ($line !~ /^\S+/) {
				$flag=2;
				#print "got here\n";#test
			}
		}
		if ($flag==2 or $flag==3) {
			if ($line =~ /^\w{2}\t\S+/) {
				$line =~ /^(\w{2})\t(\S+)/;
				$taxon = $1;
				$seq = $2;
				#print "taxon: $taxon\tseq: $seq\n";#test
				$ortho_aln{$taxon} = $seq;
				$flag=3;
				@line=();
				$taxon=();
				$seq=();
			}
		}
		if ($flag==3) {
			if ($line !~ /^\S+/) {
				while (my ($key,$value) = each (%ortho_aln) ) {
					$ntax++;
					@value = split(//,$value);
					$nchar = scalar(@value);
				}
				print OUT "$ntax\t$nchar\n";#num taxa
				$ntax=1;
				while (my($key,$value) = each(%ortho_aln)) {
					print OUT "$key        $value\n";
				}
				close OUT;
				%ortho_aln=();
				$flag=4;
				$ntax=();
				$nchar=();
				@value=();
			}
		}
		elsif ($flag==4) {
			if ($line =~ /^\[\d{8}\]/) {
				$line =~ /^\[(\d{8})\]/;
				$ortho_id = $1;

				$filename = "$ortho_id".".phy";
				open (OUT,">>",$filename) || die "Error cannot write to $filename: $!\n";
				$flag=1;
			}
		}
	}
	$i++;
}
