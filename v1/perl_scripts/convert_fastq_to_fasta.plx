#!/usr/bin/perl
#March 31, 2014 by Terri Porter
#Script to convert a directory of *.fastq.gz files into fasta format for Shadi Barcoding Illumina malaise project
#Adapt perl one-liner originally written by Robert Schmieder in the Edwards Bioinformatics Lab at San Diego State University
#USAGE perl convert_fastq_to_fasta.plx

use strict;
use warnings;

#declare var
my $i=0;
my $file;
my $fastq;
my $fasta;
my $output;
my $filesize;
my $j=0;
my $line;

#declare array
my @files;
my @in;

@files = qx(ls | grep .fastq.gz);

while ($files[$i]) {
	$file = $files[$i];
	chomp $file;

	$filesize = -s $file; #*.fastq.gz

	if ( $filesize > 100) {

		#gunzip *.fastq.gz to *.fastq
		$output = qx(gunzip $file);
		$file =~ s/\.gz//;
		$fastq = $file; #infile

		open (IN, "<", $fastq) || die "Error cannot open $fastq: $!\n";
#		@in = <IN>;
#		close IN;
		
		$fastq =~ s/\.fastq/\.fasta/;
		$fasta = $fastq; #outfile

		open (OUT, ">>", $fasta) || die "Error cannot open $fasta: $!\n";

		while (<IN>) {
			$line = $_;
			chomp $line;

			if ($line =~ /^\@/ && $j==0) {
				$line =~ s/^\@/\>/;
				print OUT $line."\n";
			}
			elsif ($j==1) {
				print OUT $line."\n";
				$j=-3;
			}
			$j++;
			$line=();
		}
		$j=0;
		close IN;
		close OUT;
	}
	
	$i++;
	$file=();
	$output=();
	$fastq=();
	$fasta=();
	@in=();
	$filesize=();

}
$i=0;





=cut

### original perl one-liner ###
### http://edwards.sdsu.edu/labsite/index.php/robert/289-how-to-convert-fastq-to-fasta ###

cat file.fastq | perl -e '$i=0; while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){print;$i=-3)$i++;}' > file.fasta

=pod
