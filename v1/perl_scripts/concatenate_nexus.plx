#!/usr/bin/perl
#October 7, 2009
#Script to concatenated reformatted gblock alignments into a single nexus file for PAUP
#USAGE: $perl concatenate_nexus.plx /path/to/infiles

use strict;
use warnings;

#declare some variables
my $dir;
my $line;
my $ortho;
my $i=0;#keep track of open files
my $j=0;
my $infile;
my $taxon;
my $seq;
my $repseq;
my $nchar_ortho;
my $nchar_cumulative=0;
my $range_start=1;

#declare some arrays
my @reformatted_gblock_files;
my @orthos;
my @seqs;
my @nchar_ortho;

#read a directory of file names into an array
$dir = "$ARGV[0]";
opendir DH, $dir;
@reformatted_gblock_files = readdir(DH);
close DH;

#at unix prompt, do $ls -a to see if the . and .. files are present, if so then:
#pop(@reformatted_gblock_files);#to remove . from end of array
#pop(@reformatted_gblock_files);#to remove ..

#grab orthologous IDs
foreach (@reformatted_gblock_files){
	$line = $_;
	if (/^\d{8}/){
		$line =~ /^(\d{8})\.+/;
		$ortho = $1;
		push (@orthos, $ortho);#add to end of array
		$i++;
	}
}
$i=2;#reset counter-skip over files . and ..
$j=0;
while ($reformatted_gblock_files[$i]){

	$infile = $ARGV[0].$reformatted_gblock_files[$i];#add path to filename
	open (IN,"<",$infile) || die "Cannot open $infile: $!";
	@seqs= ();#make sure it's empty before adding new seqs
	while (<IN>){
		$line = $_;
		chomp($line);
		if (/^\S{2}\t/){
			$line =~ /^(\S{2})\t(\S+)/g;
			$taxon = $1;
			$seq = $2;
			push (@seqs, $seq);#add to end of array
		}
	}

	$repseq = pop(@seqs);#select a representative sequence from the block (from end of array)
		@nchar_ortho = split(//,$repseq);
		$nchar_ortho = scalar(@nchar_ortho);

	#create a temp file to hold ortho ID and nchar
	open (OUT,">>","ortho_nchar\.tmp") || die "Can't write to ortho_nchar.tmp: $!";

	print OUT $orthos[$j]."\t".$nchar_ortho."\n";

	#@seqs=();#clear array
	@nchar_ortho=();
	close OUT;
	close IN;
	$i++;
	$j++;
}
$i=2;#reset counter-skip over files . and ..
$j=0;
#get nchar_cumulative
open (IN2,"<","ortho_nchar\.tmp") || die "Can't read from ortho_nchar.tmp: $!";

#print another temp file to hold ortho ID, nchar, cumulative nchar
open (OUT2,">>","ortho_nchar_cum\.tmp") || die "Can't write to ortho_nchar_cum.tmp: $!";

while (<IN2>){
	$line = $_;
	chomp($line);
	if (/\d+\t\d+/){
		$line =~ /(\d+)\t(\d+)/;
		$ortho = $1;
		$nchar_ortho = $2;
	}
	$nchar_cumulative = $nchar_cumulative + $nchar_ortho;
	print OUT2 $ortho."\t".$nchar_ortho."\t".$nchar_cumulative."\n";	
}

close IN2;
close OUT2;

#remove tempfile
unlink("ortho_nchar\.tmp");

#print a nexus formatted gblocks file
open (OUT3,">>","gblocks\.nex") || die "Can't write to gblocks.nex: $!";

print OUT3 "#NEXUS\n[Written by concatenate_nexus.plx]\n\n";
print OUT3 "Begin data;\n";
print OUT3 "Dimensions ntax=XX nchar=$nchar_cumulative;\n";#change ntax manually if necessary
print OUT3 "Format datatype=protein interleave gap=-;\n";
print OUT3 "Matrix\n\n";

while ($reformatted_gblock_files[$i]){
        $infile = $ARGV[0].$reformatted_gblock_files[$i];#add path to filename
	open (IN3,"<",$infile) || die "Cannot open $infile: $!";
        
	print OUT3 "[".$orthos[$j]."]"."\n\n";
		
	while (<IN3>){
		print OUT3;#print reformatted gblock file to nexus file
	}
	print OUT3 "\n";
	close IN3;
	$i++;
	$j++;
}

print OUT3 ";\nend;\n\n";
print OUT3 "Begin assumptions;\n";

#reopen ortho_nchar_cum.tmp file
open (IN4,"<","ortho_nchar_cum\.tmp") || die "Can't read from ortho_nchar_cum.tmp: $!";

while (<IN4>){
	$line = $_;
	chomp($line);
	if (/\d+\t\d+\t\d+/){
		$line =~ /(\d+)\t(\d+)\t(\d+)/;
		$ortho = $1;
		$nchar_ortho = $2;
		$nchar_cumulative = $3;
	}
	print OUT3 "charset ".$ortho." = ".$range_start."-".$nchar_cumulative.";\n";
	$range_start = $nchar_cumulative + 1;
 }
print OUT3 "end;\n";
close IN4;
close OUT3;

#clean up by removing old tempfiles
unlink("ortho_nchar_cum\.tmp");
		


