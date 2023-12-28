#!/usr/bin/perl
#April 11, 2014 by Terri Porter
#Script to sort R1 and R2 by reverse-complemented barcodes in the I1 .fastq.gz
#usage perl sort_R1_R2_by_I1.plx rc_barcodes.txt I1.fastq.gz R1.fastq.gz R2.fastq.gz

use strict;
use warnings;

#declare var
my $indexFile;
my $i=0;
my $line;
my $seqline;
my $sample;
my $barcode;
my $id;
my $original;
my $new;
my $prefix;
my $scalar;

#declare array
my @barcodes; #reverse complemented! tab delimited: sample\trc_barcode\n
my @index; #fastq.gz
my @r1; #fastq.gz
my @r2; #fastq.gz
my @line;
my @output;

#declare hash
my %r1; #key = idline, value=sequence | phred lines
my %r2;
my %barcode; #key=barcode, value=hash
my %sample; #key=sample, value= pipe delimited string of Illumina id lines

#collect arguments
print "\nReading infiles...\t";
open (BARCODE, "<", $ARGV[0]) || die "Error cannot open barcode file: $!\n";
@barcodes = <BARCODE>;
$scalar = scalar(@barcodes);
print "Found $scalar barcodes...\t";
close BARCODE;

$indexFile = $ARGV[1];
chomp $indexFile;

open (R1, "zcat $ARGV[2] |") || die "Error cannot open R1 file: $!\n";
@r1 = <R1>;
$scalar = scalar(@r1);
print "R1 file has $scalar lines...\t";
close R1;

open (R2, "zcat $ARGV[3] |") || die "Error cannot open R2 file: $!\n";
@r2 = <R2>;
$scalar = scalar(@r2);
print "R2 file has $scalar lines...\t";
close R2;
print "done.\n";

#hash R1 and R2 files
print "\nHashing read files...\t";
%r1 = hash_read_file(\@r1);
$scalar = keys %r1;
print "R1 hash has $scalar keys (id lines)...\t";
%r2 = hash_read_file(\@r2);
$scalar = keys %r2;
print "R2 hash has $scalar keys (id lines)...\t";
print "done.\n";

#hash barcodes
print "\nHashing barcode file...\t";
$i=0;
while ($barcodes[$i]) {
	$line = $barcodes[$i];
	chomp $line;
	
	@line = split(/\t/,$line);
	$sample = $line[0];
	$barcode = $line[1];
	$barcode{$barcode}=$sample;

	$i++;
	$line=();
	@line=();
	$sample=();
	$barcode=();
}
$i=0;
$scalar = keys %barcode;
print "Barcode hash has $scalar keys (barcodes)...\t";
print "done.\n";

#grep illumina ids from index file for each barcode
print "\nParsing illumina ids from index file...\t";
while ( ($barcode,$sample) = each (%barcode) ) {

	@output = qx(zcat $indexFile | grep $barcode -B 1);
	
	$i=0;
	while ($output[$i]) {
		$line = $output[$i];
		chomp $line;

		if ($line =~ /^@/) {
#			$line =~ s/^@//;
			$id = $line; #illumina fastq header including @
			if (exists $sample{$sample}) {
				$original = $sample{$sample};
				$new = $original."|".$id;
				$sample{$sample} = $new;
			}
			else {
				$sample{$sample} = $id;
			}
		}
		else {
			$i+=4; #increment through file faster
			next;
		}

		$i++;
		$line=();
		$id=();
		$original=();
		$new=();
	}
	$i=0;

	@output=();

}
print "Sample hash has $scalar keys (samples)...\t";
print "done.\n";

#parse through R1/R2 fastq.gz files and create fastq files for each sample
print "\nParsing R1 by sample...\t";
$prefix = "R1";
parse_fastq(\%r1,\$prefix);
print "done.\n";

print "\nParsing R2 by sample...\t";
$prefix = "R2";
parse_fastq(\%r2,\$prefix);
print "done.\n";

####################

sub hash_read_file {

#declare var
my $line;
my $j;
my $seqline;
my $k;
my $phredline;

#declare array
my @infile = @{$_[0]}; #fastq

#declare hash
my %hash;

while ($infile[$i]) {
	$line = $infile[$i];
	chomp $line;

	if ($line =~ /^@/) {
		$id = $line;
		$j=$i+1; #get seq line
		$seqline = $infile[$j];
		chomp $seqline;
		$k = $i+3; #get phred line
		$phredline = $infile[$k];
		chomp $phredline;
		$hash{$id} = $seqline."|".$phredline;
		$i=$k;
	}
	$i++;
	$line=();
	$id=();
	$j=();
	$seqline=();
	$k=();
	$phredline=();

}
$i=0;

return %hash;

}


####################

sub parse_fastq {

#declare var
my $idline;
my $i=0;
my $id;
my $seqline;
my $phredline;
my $prefix = ${$_[1]}; #R1 or R2
my $outfile;
my $value;

#declare array
my @idline;
my @value;

#declare hash
my %readfile = %{$_[0]}; #key = id, value = seqline | phredline

print "\nGot to parse_fastq()\n";

while ( ($sample,$idline) = each (%sample)) {
	print "Trying to process sample: $sample\tidline: $idline\n";
	@idline = split(/\|/,$idline);

	$outfile = $sample.".".$prefix.".fastq";
	open (OUT, ">>", $outfile) || die "Error cannot open $outfile\n";
	
	$i=0;
	while ($idline[$i]) {
		$id = $idline[$i];

		if ($prefix eq 'R2') { 
			$id =~ s/1:N:0:0/2:N:0:0/;
		}

		if (exists $readfile{$id}) {
			$value = $readfile{$id};
			@value = split(/\|/,$value);
			$seqline = $value[0];
			$phredline = $value[1];
			print OUT "$id\n$seqline\n+\n$phredline\n";
		}
		else {
			print "Cannot find id $id for sample $sample in $prefix file\n";
		}
		
		$i++;
		$id=();
		$value=();
		@value=();
		$seqline=();
		$phredline=();

	}
	$i=0;
	close OUT;
	@idline=();
	$outfile=();

}


}

####################
