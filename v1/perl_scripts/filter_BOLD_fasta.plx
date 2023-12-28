#!/usr/bin/perl
#Apr. 26/18 script to filter fasta from BOLD Public Data API using query_BOLD.sh taxonlist.txt
#Remove non-COI seqs, don't worry about GB-SUPPRESSED
#USAGE perl filter_BOLD_fasta.plx cat.fasta

#declare var
my $outfile="cat.fasta.filtered";
my $line;
my $i=0;
my $gene;

#declare array
my @fasta;
my @line;

#read in fasta file
open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@fasta = <IN>;
close IN;

#create filtered fasta outfile
open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

#parse through the file and filter out unwanted records
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	#fasta header found
	if ($line =~ /^>/) {
		@line = split(/\|/,$line); #[0] BOLD ID, [1] Taxon name, [2] gene, [3] GB accession
		$gene = $line[2];

		#check for non-COI gene
		if ($gene !~ /COI/) {
			$i+=2; #skip two lines to go to the next fast header
		}
		else {
			print OUT $line."\n";
		}
	}
	else {
		print OUT $line."\n";
	}
	$i++;
	$line=();
	@line=();
	$gene=();
}
$i=0;
