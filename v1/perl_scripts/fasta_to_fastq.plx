#!/usr/bin/perl
#March 5, 2014 by Terri Porter
#Script to turn a fasta file into a fastq file, to enable FC and BR pairing with SeqPrep (after all normal processing, clustering, filtering, etc. has already been done)
#usage perl fasta_to_fastq.plx FC.fasta BR.fasta

use strict;
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;

#declare var
my $i=0;
my $j;
my $line;
my $header;
my $seq;
my $length;
my $nextline;
my $rc; #revese-complement
my $seqin;
my $seqobj;
my $seqout;
my $flag=0;

#declare array
my @rc;
my @seq;
my @temp;


#process FC or forward reads first
open (OUT1, ">>", "fwd.fastq") || die "Error cannot open fwd.fastq: $!\n";

$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$ARGV[0]);

while ($seqobj = $seqin->next_seq()) {

	$header = $seqobj->display_id;
	$seq = $seqobj->seq;
	@seq = split(//,$seq);
	$length = scalar(@seq);

	print OUT1 "\@$header\n"; #tested, SeqPrep does not require illumina header format, just a unique header
	print OUT1 "$seq\n";
	print OUT1 "+\n";
	print OUT1 'I' x $length;
	print OUT1 "\n";

	$header=(); 
	$seq=();
	@seq=();
	$length=();
}
close OUT1;
$seqin=();

#process BR or reverse reads second, don't forget to reverse-complement!
$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$ARGV[1]);
$seqout = Bio::SeqIO->new(-format=>'Fasta', -file=>'>>temp.fa');

while ($seqobj = $seqin->next_seq()) {

	$header = $seqobj->display_id;
	$rc=$seqobj->revcom;
	$seqout->write_seq($rc);

	$header=();
	$rc=();
	$seqout=();
}
$seqin=();
$seqout=();

#parse the temp.fa file to change header format and remove newline breaks in seq part
open (TEMP, "<", 'temp.fa') || die "Error cannot open temp.fa: $!\n";
@temp = <TEMP>;
close TEMP;

open (OUT2, ">>", "rev.fastq") || die "Error cannot open rev.fastq: $!\n";

while ($temp[$i]) {
	$line = $temp[$i];
	chomp $line;

	if ( $flag==0 && $line =~ /^>/) {
		$header = $line;
		$header =~ s/>//;
		$flag=1;	
		print OUT2 "\@$header\n";
	}
	elsif ( $flag==1 && $line !~ /^>/ ) {
		$seq = $line;
		$flag=2;
	}
	elsif ($flag==2 && $line !~ /^>/) {
		$seq = $seq.$line;
		$flag=2;
	}
	elsif ($flag==2 && $line =~ /^>/) {	
		print OUT2 "$seq\n";
		print OUT2 "+\n";
		@seq = split(//,$seq);
		$length = scalar(@seq);
		print OUT2 'I' x $length;
		print OUT2 "\n";
		$i--;
		$flag=0;
		@seq=();
		$seq=();
		$length=();
	}
	$i++;

	$line=();
	$header=();
}

#don't forget to process the last sequence!
print OUT2 "$seq\n";
print OUT2 "+\n";
@seq = split(//,$seq);
$length=scalar(@seq);
print OUT2 'I' x $length;
print OUT2 "\n";
$flag=0;
@seq=();
$seq=();
$length=();
$i=0;
close OUT2;
