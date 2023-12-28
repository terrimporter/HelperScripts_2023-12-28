#!/usr/bin/perl
#April 21, 2014 by Terri Porter
#Script to run usearch on ITS2 extracted files (Nagissa's Caspian and Mediterranean sea samples)
#USAGE perl Nagissa_run_usearch.plx

use strict;
use warnings;
use Cwd;
use File::Copy;

#declare var
my $i=0;
my $file;
my $file2;
my $jobs = 20; #number of blast files to run at the same time

#my $original_dir; #where the ITS2 extracted fastas are
#my $NBCpath = "/home/terri/rdp_classifier2.6";
#my $cwd;

#declare array
my @output;

@output = qx(ls | grep fasta);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	$file2 = run_usearch(\$file);

#	$original_dir = cwd();
#	print "\n$cwd\n";

#	chdir $NBCpath;

#	$cwd = cwd();
#	print "$cwd\n";

	run_BLAST(\$file2);

	$file=();
	$file2=();
	$i++;
}
$i=0;

####################

sub run_usearch {

#declare var
my $infile = ${$_[0]}; #ITS2 extracted fasta
my $result;
my $prefix;
my $filename;
my $derep;
my $sorted;
my $nonch_outfile;
my $ch_outfile;
my $sorted2;
my $centroids;
my $centroids2;

#declare array
my @infile;

@infile = split(/\./,$infile);
$filename = $infile[0];
$filename =~ s/_paired//;
$prefix = $filename;

$derep = $prefix.".derep.fa";

$result = system("usearch6.0.307_i86linux32 -derep_fulllength $infile -sizeout -output $derep");

$sorted = $prefix.".sorted.fa";

$result = system("usearch6.0.307_i86linux32 -sortbysize $derep -output $sorted");

$nonch_outfile = $prefix.".nonch.fa";
$ch_outfile = $prefix.".ch.fa";

$result = system("usearch6.0.307_i86linux32 -uchime_denovo $sorted -nonchimeras $nonch_outfile -chimeras $ch_outfile");

$sorted2 = $prefix.".sorted2.fa";

$result = system("usearch6.0.307_i86linux32 -sortbysize $nonch_outfile -output $sorted2");

$centroids = $prefix.".centroids.fa";

$result = system("usearch6.0.307_i86linux32 -cluster_smallmem $sorted2 -id 0.97 -sizein -sizeout -centroids $centroids");

$centroids2 = $prefix.".centroids_minsize2.fa";

$result = system("usearch6.0.307_i86linux32 -sortbysize $centroids -minsize 2 -output $centroids2");

return $centroids2;

}

####################

sub run_BLAST {

#declare var
my $file = ${$_[0]}; #minsize 2 file from usearch
my $outfile = $file.".blastn";
#my $result;
my $prefix;
my $original_dir;
my $output;
my $old_file;
my $new_file;
my $cwd;
my $file2;

#declare array
my @file;

@file = split(/\./,$file);
$prefix = $file[0];

$original_dir = cwd();
print "\nCWD: $original_dir\n";
mkdir $prefix;
chdir $prefix;
$cwd = cwd();
print "CWD: $cwd\n";

$old_file = $original_dir."/".$file;
$new_file = $original_dir."/".$prefix."/".$file;

copy $old_file, $new_file;

$file2 = reformat_fasta(\$file);

$output = qx(split -l 1000 $file2);

$output = qx(ls | grep '^x' | parallel -j $jobs "blastn -task megablast -db /1/scratch/blastdb/nt -query {} -out {}.blastn -evalue '1e-20' -outfmt 0 -num_descriptions 100 -num_alignments 100");

chdir $original_dir;
$cwd = cwd();
print "CWD: $cwd\n";

}

####################

sub reformat_fasta {

#declare var
my $file = ${$_[0]};
my $prefix;
my $j=0;
my $line;
my $header;
my $flag=0;
my $k;
my $scalar;
my $outfile;
my $seq;
my $nextline;
my $maxline;
my $seq2;

#declare array
my @in;
my @file;

open (IN, "<", $file) || die "Error cannot open $file\n";
@in = <IN>;
close IN;

@file = split(/\./,$file);
$prefix = $file[0];
$outfile = $prefix.".reformatted";

open (OUT, ">>", $outfile) || die "Error cannot open $outfile:$!\n";

while ($in[$j]) {
	$line = $in[$j]; #do not chomp!
	
	if ($j==0 && $line =~ /^>/) {
		$line =~ s/=/_/g;
		print OUT $line; #header
	}
	elsif ($j > 0 && $line !~ /^>/) {
		chomp $line; #seq
		print OUT $line;
	}
	elsif ($j > 0 && $line =~ /^>/) {
		$line =~ s/=/_/g;
		print OUT "\n";
		print OUT $line; #header
	}
	$j++;

}
$j=0;
$line=();
close OUT;

return $outfile;

}

####################
