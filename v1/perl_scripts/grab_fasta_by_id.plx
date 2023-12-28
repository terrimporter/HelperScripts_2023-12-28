#!/usr/bin/perl
#Oct.3, 2013 by Terri Porter
#Script to grab the associated fasta (from benthos_pipeline.plx) for each readid assigned by NBC (from parse_NBC.plx)
#USAGE perl grab_fasta_by_id.plx NBC.parsed cluster_98otus_minsize2.fasta.newlineremoved

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $id; #readid including size from usearch...
my $readid;
my $size;
my $kingdom;
my $kingdom_bp;
my $phylum;
my $phylum_bp;
my $class;
my $class_bp;
my $order;
my $order_bp;
my $family;
my $family_bp;
my $genus;
my $genus_bp;
my $lineage;
my $j;
my $nextline;

#declare array
my @readid;
my @id;
my @fasta;
my @line;

#declare hash
my %header; #indexed by readid

open (READID, "<", $ARGV[0]) || die "Error cannot open NBC.parsed: $!\n";
@readid = <READID>;
close READID;

open (FASTA, "<", $ARGV[1]) || die "Error cannot open cluster_98otus_minsize2.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

#hash readids and lineage
while ($readid[$i]) {
	$line = $readid[$i];
	chomp $line;

	if (length($line) > 0) {
		@line = split(/\t/,$line);
		$id = $line[0];
#		@id = split(';',$id);
#		$readid = $id[0];
#		$readid =~ s/^>//;
#		$size = $id[1];
#		$size =~ s/size=//;
		$kingdom = $line[1];
		$kingdom_bp = $line[2];
		$phylum = $line[3];
		$phylum_bp = $line[4];
		$class = $line[5];
		$class_bp = $line[6];
		$order = $line[7];
		$order_bp = $line[8];
		$family = $line[9];
		$family_bp = $line[10];
		$genus = $line[11];
		$genus_bp = $line[12];
		$lineage = $kingdom."|".$kingdom_bp."|".$phylum."|".$phylum_bp."|".$class."|".$class_bp."|".$order."|".$order_bp."|".$family."|".$family_bp."|".$genus."|".$genus_bp;
#		$header{$readid} = $readid."|".$size."|".$lineage."\n";
		$header{$id} = $lineage;
	}
	$i++;
	$line=();
	@line=();
	$id=();
#	$readid=();
#	$size=();
	$kingdom=();
	$kingdom_bp=();
	$phylum=();
	$phylum_bp=();
	$class=();
	$class_bp=();
	$order=();
	$order_bp=();
	$family=();
	$family_bp=();
	$genus=();
	$genus_bp=();
	$lineage=();
}
$i=0;

#search fasta for just the ones with a good taxonomic assignment

open (OUT, ">>", "NBC.parsed.fasta") || die "Error cannot open NBC.parsed.fasta:$!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$id = $1;

		if (exists $header{$id}) {
			@id = split(';',$id);
			$readid = $id[0];
			$readid =~ s/^>//;
			$size = $id[1];
			$size =~ s/size=//;
			$lineage = $header{$id};
			print OUT ">$readid|$size|$lineage\n";
			$j=$i+1;
			$nextline = $fasta[$j];
			chomp $nextline;
			print OUT "$nextline\n";
			$i+=2;
			$id=();
			@id=();
			$readid=();
			$size=();
			$lineage=();
			$j=();
			$nextline=();
			next;
		}
		else {
			print "$id does not have a good taxonomic assignment, so not printed.\n";
		}
	}
	$i++;
	$line=();
}
$i=0;

close OUT;
