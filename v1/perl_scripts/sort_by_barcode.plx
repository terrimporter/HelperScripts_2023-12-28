#!/usr/bin/perl
#March 25, 2014 by Terri Porter
#Script to process Nagissa's Knoxville ITS Illumina reads from Caspian and Mediterranean Seas
#usage perl sort_by_barcoede.plx file.fasta file.qual sample_barcode.txt
#First, be sure to convert file.fastq.gz into a file.fasta for simpler manipulation
#sample_barcode.txt is a tab delimited file: sample\tbarcode
#The sample name, first five characters should be unique, no spaces, underscores ok

use strict;
use warnings;
use Cwd;

#declare var
my $sample;
my $barcode;

#declare array
my @fasta;
my @qual;
my @map;
my @simpleFasta;
my @simpleQual;

#declare hash
my %barcode_sample; #key = barcode, value = sample (5 letters)
my %sample_readid; #key = sample, value = readid (5 letters)
my %illuminareadid_readid; #key = illuminareadid, value = readid

get_infiles(); #fasta, qual, barcode_sample.map

setup_folders();

simplify_readids(); #in .fasta and .qual, also create illuminareadid_readid.map

sort_fasta_by_barcode(); #create sample_readid.map

create_qual_by_fasta();

#run_usearch(); #dereplicate with sizeout, sort by cluster size, cluster with sizein sizeout, remove singletons

#create_maps(); #Readid-derep.map, Readid-cluster.map

#print "\nNow ready to manually run FungalITSextractor\n";

####################

sub get_infiles {

open (FASTA, "<", $ARGV[0]) || die "Error can't open fasta file: $!\n";
@fasta = <FASTA>;
close FASTA;

open (QUAL, "<", $ARGV[1]) || die "Error can't open qual file: $!\n";
@qual = <QUAL>;
close QUAL;

open (MAP, "<", $ARGV[2]) || die "Error can't open sample_barcode.map: $!\n";
@map = <MAP>;
close MAP;

}

####################

sub setup_folders {

#declare var
my $i=0;
my $line;
my $directory;

#declare array
my @line;

#parse sample_barcode.map to get folder names (first five letters of sample column)
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$sample = $line[0];
	$barcode = $line[1];
	$barcode_sample{$barcode} = $sample;

	$i++;
	$line=();
	@line=();
	$sample=();
	$barcode=();
}
$i=0;

#create a new directory for each sample
while ( ($barcode,$sample) = each (%barcode_sample) ) {
	$directory = $sample;
	
	unless (-e $directory or mkdir $directory) {
		die "Unable to create $directory\n";
	}
	$directory=();
}

}

####################

sub simplify_readids {

#declare var
my $i=0;
my $illuminaReadid;
my $simple_readid;
my $readid_counter = 1;
my $newline;
my $line;

open (MAP2, ">>", "readid_illuminareadid.map") || die "Error cannot open readid_illuminareadid.map: $!\n";

open (OUT, ">>", "simple.fasta") || die "Error cannot open simple.fasta: $!\n";

#first, process the fasta file ### assumes no line-wrapping within seq line
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$illuminaReadid = $1;
		$simple_readid = $readid_counter;

		print MAP2 "$simple_readid\t$illuminaReadid\n";
		$illuminareadid_readid{$illuminaReadid} = $simple_readid;

		$readid_counter++;
		$newline = ">$simple_readid";
		print OUT "$newline\n";
	}
	else {
		print OUT "$line\n";
	}
	$i++;
	$line=();
	$illuminaReadid=();
	$simple_readid=();
	$newline=();

}
$i=0;
close OUT;

#second, process the qual file ### assumes no line wrapping within phredline
$i=0;

open (OUT2, ">>", "simple.qual") || die "Error cannot open simple.qual: $!\n";

while($qual[$i]) {
	$line = $qual[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$illuminaReadid = $1;

		if (exists $illuminareadid_readid{$illuminaReadid}) {
			$simple_readid = $illuminareadid_readid{$illuminaReadid};

			$newline = ">$simple_readid";
			print OUT2 "$newline\n";
		}
	}
	else {
			print OUT2 "$line\n";
	}

	$i++;
	$line=();
	$illuminaReadid=();
	$simple_readid=();
	$newline=();

}
$i=0;
close OUT2;

}

####################

sub sort_fasta_by_barcode {

#declare var
my $i=0; 
my $line;
my $readid;
my $seq;
my $sub; #substring  = should be barcode
my $original;
my $new;
my $readidline;
my $newfile;
my $output;

#declare array
my @readidline;
my @output;

open (FASTA, "<", "simple.fasta") || die "Error can't open simple.fasta: $!\n";
@simpleFasta = <FASTA>;
close FASTA;

#parse fasta to create readid_sample.map
while($simpleFasta[$i]) {
	$line = $simpleFasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\S+)/;
		$readid = $1;
	}
	else { ### assumes no line-wrapping in seq
		$seq = $line;
		$sub = substr $seq, 0, 12; #grab first 12 bp of seq, this should be the barcode
		if (exists $barcode_sample{$sub}) {
			$sample = $barcode_sample{$sub};

			if (exists $sample_readid{$sample}) {
				$original = $sample_readid{$sample};
				$new = $original."|".$readid;
				$sample_readid{$sample} = $new;
			}
			else {
				$sample_readid{$sample}=$readid;
			}
		}
	}
	
	$i++;
	$line=();
#	$readid=();
	$seq=();
	$sub=();
	$sample=();

}
$i=0;

#for each sample, print appropriate readids to new fasta file
while ( ($barcode,$sample) = each (%barcode_sample) ) {
	$newfile = $sample."/".$sample.".fasta"; #path to sample/
	open (OUT, ">>", $newfile) || die "Error cannot open $newfile:$!\n";

	if (exists $sample_readid{$sample}) {
		$readidline = $sample_readid{$sample};
		@readidline = split(/\|/,$readidline);

		while ($readidline[$i]) {
			$readid = $readidline[$i];

			@output = qx(grep '\>$readid\$\' simple.fasta -A 1); # grep -A prints lines after match too

			foreach $output (@output) {
				chomp $output;
				print OUT $output."\n";
#				print $output."\n";#test
			}

			$i++;
			$readid=();
			@output=();
			$output=();
		}
		$i=0;
		close OUT;

	}
	$readidline=();
	@readidline=();

}

}

####################

sub create_qual_by_fasta {

#declare var
my $fastafile;
my $i=0;
my $line;
my $readid;
my $qualfile;
my $output;

#declare array
my @in;
my @output;

while ( ($barcode,$sample) = each (%barcode_sample) ) {
	$fastafile = $sample."/".$sample.".fasta";

	if (-e $fastafile) {
		open (FASTA, "<", $fastafile) || die "Error cannot open $fastafile: $!\n";
		@in = <FASTA>;
		close FASTA;

		$qualfile = $sample."/".$sample.".qual";
		open (QUAL, ">>", $qualfile) || die "Error cannot open $qualfile: $!\n";

		while ($in[$i]) {
			$line = $in[$i];
			chomp $line;

			if ($line =~ /^>/) {
				$line =~ /^>(\S+)/;
				$readid = $1;

				@output = qx(grep '\>$readid\$\' -A 1 simple.qual);

				foreach $output (@output) {
					chomp $output;
					print QUAL "$output\n";
				}
				$readid=();
				@output=();
				$output=();
			}

			$i++;
			$line=();
		}
		$i=0;
		close QUAL;
	}
}

}

####################

sub run_usearch {

#declare var
my $originaldir;
my $cwd;
my $output;
my $infile;
my $outfile;
my $outfile_uc;
my $minseqlength = 80; ##### change minimum seq length here #####
my $sorted;
my $nonch_outfile;
my $ch_outfile;
my $sorted2;
my $centroids;
my $centroids_uc;
my $minsize;


$originaldir = cwd;
print "originaldir: $originaldir\n"; #test

while ( ($barcode,$sample) = each (%barcode_sample) ) {

	chdir $sample;
	$cwd = cwd;
	print "cwd: $cwd\n"; #test

	$infile = $sample.".fasta";
	$outfile = $infile.".derep";
	$outfile_uc = $infile.".derep.uc";
	#currently using usearch v7.0.959_i86linux32
	$output = system("usearch -minseqlength $minseqlength -derep_prefix $infile -sizeout -output $outfile -uc $outfile_uc");

	$sorted = $outfile.".sorted";
	$output = system("usearch -sortbysize $outfile -output $sorted");

	$nonch_outfile = $sorted.".nonch";
	$ch_outfile = $sorted.".ch";
	$output = system("usearch -uchime_denovo $sorted -nonchimeras $nonch_outfile -chimeras $ch_outfile");

	$sorted2 = $nonch_outfile.".sorted";
	$output = system("usearch -sortbysize $nonch_outfile -output $sorted2");

	$centroids = $sorted2.".centroids";
	$centroids_uc = $sorted2.".centroids.uc";
	$output = system("usearch -cluster_smallmem $sorted2 -id 0.97 -centroids $centroids -uc $centroids_uc -sizein -sizeout -usersort");

	$minsize = $centroids.".minsize2";
	$output = system("usearch -sortbysize $centroids -minsize 2 -output $minsize");

	$cwd=();
	$infile=();
	$outfile=();
	$outfile_uc=();
	$output=();
	$sorted=();
	$nonch_outfile=();
	$ch_outfile=();
	$sorted2=();
	$centroids=();
	$centroids_uc=();
	$minsize=();


	chdir $originaldir;
	$cwd = cwd;
	print "cwd: $cwd\n"; #test

}

}

####################
