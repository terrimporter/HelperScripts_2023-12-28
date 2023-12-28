#!/usr/bin/perl
#Feb. 9, 2012 by Terri Porter
#Script to use a set of OTU ids from MEGAN, ex. Inocybe maculata, to grab all associated read IDS from results.clstr, then grab all sequences from All_longer.fasta
#usage perl grab_subset.plx OTUid.list results.clstr All_longer.fasta

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $MEGAN_OTU;
my $readIDline;
my $CLSTR_OTU;
my $CLSTR_readid;
my $readid;
my $j;
my $next_line;
my $seq;

#declare array
my @OTU;
my @CLSTR;
my @FASTA;
my @line;
my @MEGAN_OTUs;

#declare hash
my %readClusterMap1;
my %readClusterMap2;
my %fasta;

open (OTU,"<", $ARGV[0]) || die "Error cannot read OTU list file: $!\n";
@OTU = <OTU>;
close OTU;

open (CLUSTER,"<",$ARGV[1]) || die "Error cannot read results.clstr file: $!\n";
@CLSTR = <CLUSTER>;
close CLUSTER;

open (FASTA,"<",$ARGV[2]) || die "Error cannot read fasta file: $!\n";
@FASTA = <FASTA>;
close FASTA;

#grab OTUs from MEGAN tab-delimited file

while ($OTU[$i]) {
	$line = $OTU[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$MEGAN_OTU = $line[0];
	push (@MEGAN_OTUs, $MEGAN_OTU);

	$i++;
	$line=();
	@line=();
	$MEGAN_OTU=();

}
$i=0;

#print "MEGAN OTU list:@MEGAN_OTUs\n";

#for each OTU in results.clstr, grab all associated reads into an OTU_readid map

while ($CLSTR[$i]) {
	$line = $CLSTR[$i];
	chomp $line;

	if ($line =~ /\*/) {
		@line = split(/\s+/,$line);

		$readIDline = $line[2];
		$readIDline =~ s/^>//;
		$readIDline =~ s/\.\.\.$//;
		$CLSTR_OTU = $readIDline;
#		print "CLSTR_OTU:$CLSTR_OTU\n";

		$readIDline=();
		@line=();
	
	}

	elsif ($line =~ /%$/) {
		@line = split(/\s+/,$line);
		$readIDline = $line[2];
#		print "readIDline:$readIDline\n";
		$readIDline =~ s/^>//;
		$readIDline =~ s/\.\.\.$//;
		$CLSTR_readid = $readIDline;

		$readClusterMap1{$CLSTR_readid} = $CLSTR_OTU; #contains all OTUs and associated reads

		@line=();
		$readIDline=();
		$CLSTR_readid = ();
	}

	elsif ($line =~ /^>/) {
		$CLSTR_OTU=();
	}
	$i++;
}
$i=0;

#for each MEGAN OTU grab all associated CLSTR readids

while ($MEGAN_OTUs[$i]) {
	$MEGAN_OTU = $MEGAN_OTUs[$i];

	while (my($key,$value) = each (%readClusterMap1)) {
		if ($value eq $MEGAN_OTU) {
			$readClusterMap2{$key} = $MEGAN_OTU; #contains only MEGAN OTUs and associated reads
		}
	}
	$i++;
	$MEGAN_OTU=();
}
$i=0;

#read fasta file into a hash for faster searching

while($FASTA[$i]) {
	$line = $FASTA[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w+)/;
		$readid = $1;
		$j=$i+1;
		
		$next_line = $FASTA[$j];
		chomp $next_line;
		
		$seq = $next_line;
		$fasta{$readid} = $seq;
	}

	$i++;
	$readid=();
	$next_line=();
	$seq=();

}
$i=0;
$j=0;

#get the sequence for each readid associated with each MEGAN OTU from fasta hash

open (OUT,">>","read.fasta") || die "Error cannot write to read.fasta: $!\n";

while (my($key,$value) = each(%readClusterMap2)) {
	$readid = $key;
	$seq = $fasta{$readid};
	print OUT ">$readid\n$seq\n";
}

close OUT;
