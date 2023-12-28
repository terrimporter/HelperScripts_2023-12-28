#!/usr/bin/perl
#Terri Porter, May 29, 2013
#Script to parse fasta.trimmed files from part1_Shadi_052713.plx
#manually combine fwd and rev files for markers to be combined, ex. 16Sv3, 16Sv6
#run this individually for each marker

use strict;
use warnings;

#global var
my $infile;
my $outfile;
my $nodes;

print "Create readid_sample map then concatenate into a single file.\n";
create_map_then_cat();

print "Get stats.\n";
get_stats();

print "Run usearch to dereplicate, chimera check (both algorithms), cluster, and remove singletons.\n";
run_usearch();

print "Making centroid_readid mapping file for dereplicated reads.\n";
$infile = 'cat_derep.uc';
$outfile = 'centroid_readid_derep.map';
parse_uc(\$infile,\$outfile);

print "Making centroid_clustering mapping file.\n";
$infile = 'cat_97otus.uc';
$outfile = 'centroid_readid_OTU.map';
parse_uc(\$infile,\$outfile);

print "Creating readid_OTUcentroid mapping file.\n";
make_map();

print "Remove new lines from sequence part of cat_97otus_minsize2.fa and fix header.\n";
$infile = "cat_97otus_minsize2.fa";
$outfile = "OTUs.fa";
fix_fasta(\$infile,\$outfile);

print "Splitting OTUs.fa into smaller files then doing BLAST.\n";
$nodes = 10;
split_then_blast(\$nodes);

####################

sub fix_fasta {

#var
my $line;
my $infile = ${$_[0]};
my $outfile = ${$_[1]};
my $i=0;
my $scalar;
my $limit;

#array
my @in;

open (IN, "<", $infile) || die "Error can't open infile: $!\n";
@in = <IN>;
close IN;
$scalar = scalar(@in);
$limit = $scalar-1;

print "scalar:$scalar\n";

open (OUT,">>",$outfile) || die "Cannot open newline removed fasta file:$!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($i==0) {
		print OUT $line,"\n";
	}
	elsif ($i>1 && $i<$limit && $line !~ /^>/) {
		print OUT $line;
	}
	elsif ($i>1 && $i<$limit && $line =~ /^>/) {
		print OUT "\n",$line,"\n";
	}
	elsif ($i==$limit && $line !~ /^>/) {
		print OUT $line,"\n";
	}
	$i++;
}
$i=0;
close OUT;

}

####################

sub split_then_blast {

#var
$nodes = ${$_[0]};
my $output;

#array
my @output;

@output = qx(split -l 2000 OTUs.fa);

$output = qx(ls | grep '^x' | parallel -j $nodes "blastn -task megablast -db /1/scratch/terri/nt/nt -query {} -out {}.blastn -evalue '1e-10' -outfmt 0 -num_descriptions 100 -num_alignments 100");

}

####################

sub make_map {

#var
my $i=0;
my $line;
my $readid;
my $centroid;
my $OTUcentroid;
my $val;
my $readid2;

#array
my @map1;
my @line;
my @map2;
my @map3;

#hash
my %first; #indexed by readid
my %second; #indexed by readid
my %third; #indexed by readid2 (centroid from %second)
my %final; #indexed by readid

open (MAP1, "<", "readid_sample.map") || die "Can't open first map file: $!\n";
@map1 = <MAP1>;
close MAP1;

#initialize hash with map1 readids as keys
while ($map1[$i]) {
	$line = $map1[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$readid = $line[0];
	$first{$readid} = 1;

	$i++;
	$line=();
	@line=();
	$readid=();
}
$i=0;

open (MAP2, "<", "centroid_readid_derep.map") || die "Can't open second map file: $!\n";
@map2 = <MAP2>;
close MAP2;

#hash map 2
while ($map2[$i]) {
	$line = $map2[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$readid = $line[1];
	$readid2 = $line[0];

	$second{$readid} = $readid2;

	$i++;
	$line=();
	@line=();
	$readid=();
	$readid2=();
}
$i=0;

open (MAP3, "<", "centroid_readid_OTU.map") || die "Can't open third map file: $!\n";
@map3 = <MAP3>;
close MAP3;

#hash map3
while ($map3[$i]) {
	$line = $map3[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$readid2 = $line[1];
	$OTUcentroid = $line[0];

	$third{$readid2} = $OTUcentroid;
	
	$i++;
	$line=();
	@line=();
	$readid2=();
	$OTUcentroid=();
}
$i=0;

open (OUT, ">>", "readid_OTUcentroid.map") || die "Error can't open final map: $!\n";

while( ($readid,$val) = each(%first) ) {

	if (exists $second{$readid}) {
		$readid2 = $second{$readid};
		
		if (exists $third{$readid2}) {
			$OTUcentroid = $third{$readid2};
			print OUT "$readid\t$OTUcentroid\n";
		}
#		else {
#			print "Problem finding readid2 $readid2 from 2nd map file in 3rd map file, probably a chimera.\n";
#		}
	}
	else {
		print "Problem finding readid $readid from 1st map file in 2nd map file\n";
	}
}
close OUT;

}

####################

sub parse_uc {

#var
$infile = ${$_[0]};
$outfile = ${$_[1]};
my $i=0;
my $line;
my $readid;
my $centroid_original;
my $centroid_current;

#array
my @in;
my @line;

open (IN, "<", $infile) || die "Error can't open uc infile file: $!\n";
@in = <IN>;
close IN;

open (OUT, ">", $outfile) || die "Error can't open parsed uc map file: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;
	@line = split(/\t/,$line);

	if ($line =~ /^S/) {
		$centroid_original = $line[8];
		if ($centroid_original =~ /;/) {
			$centroid_original =~ s/;size=\d+;//;
			print OUT "$centroid_original\t$centroid_original\n";
		}
		else {
			print OUT "$centroid_original\t$centroid_original\n";
		}
#		print "centroid_original:$centroid_original\n";

	}
	elsif ($line =~ /^H/) {	
		$readid = $line[8];
		$centroid_current = $line[9];

		if ($centroid_current =~ /;/) {
			$centroid_current =~ s/;size=\d+;//;
			$readid =~ s/;size=\d+;//;
			print OUT "$centroid_current\t$readid\n";
		}
		else {
			print OUT "$centroid_current\t$readid\n";
		}
#		print "centroid_current:$centroid_current\treadid:$readid\n";
	}

	$i++;
	$line=();
	@line=();
	$readid=();
	$centroid_current=();
	$centroid_original=();

}
$i=0;
close OUT;

}

####################

sub run_usearch {

#var
my $result;

#array
my @output;

$result = system('usearch6.0.307_i86linux32 -derep_prefix cat.fasta -sizeout -output cat_derep.fa -uc cat_derep.uc');

#$result = system('usearch6.0.307_i86linux32 -cluster_smallmem ITS1_derep.fa -id 0.99 -centroids ITS1_denoised.fa -sizein -sizeout');

$result = system('usearch6.0.307_i86linux32 -sortbysize cat_derep.fa -output cat_derep.fa.sorted');

$result = system('usearch6.0.307_i86linux32 -uchime_denovo cat_derep.fa.sorted -nonchimeras cat_nonch_denovo.fa -chimeras cat_ch_denovo.fa');

#$result = system('usearch6.0.307_i86linux32 -uchime_ref cat_derep.fa.sorted -db cat_nonch_denovo.fa -strand plus -nonchimeras cat_nonch_ref.fa -chimeras cat_ch_ref.fa');
#omit this step, because reads filtered out using denovo method are removed forever and doing this step makes it doubly stringent!

$result = system('usearch6.0.307_i86linux32 -sortbysize cat_nonch_denovo.fa -output cat_sorted.fa');

$result = system('usearch6.0.307_i86linux32 -cluster_smallmem cat_sorted.fa -id 0.97 -sizein -sizeout -centroids cat_97otus.fa -uc cat_97otus.uc');

$result = system('usearch6.0.307_i86linux32 -sortbysize cat_97otus.fa -minsize 2 -output cat_97otus_minsize2.fa');

#count number of seqs in each file and print to screen for easy recording
$result = qx(grep ">" cat.fasta | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" cat_derep.fa | wc -l);
chomp $result;
print $result."\t";

#$result = qx(grep ">" ITS1_denoised.fa | wc -l);
#chomp $result;
#print $result."\t";

$result = qx(grep ">" cat_nonch_denovo.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" cat_ch_denovo.fa | wc -l);
chomp $result;
print $result."\t";

#$result = qx(grep ">" cat_nonch_ref.fa | wc -l);
#chomp $result;
#print $result."\t";

#$result = qx(grep ">" cat_ch_ref.fa | wc -l);
#chomp $result;
#print $result."\t";

$result = qx(grep ">" cat_97otus.fa | wc -l);
chomp $result;
print $result."\t";

$result = qx(grep ">" cat_97otus_minsize2.fa | wc -l);
chomp $result;
print $result."\n";

}

####################

sub get_stats {

#declare var
my $i=0;
my $file;
my $sample;
my $j=0;
my $line;
my $count=0;

#declare array
my @output;
my @in;
my @file;

@output = qx(ls | grep 'fasta.trimmed\$\');

open (OUT, ">>", "stats.txt") || die "Error cannot open outfile: $!\n";

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	open (IN, "<", $file) || die "Error cannot open $file: $!\n";
	@in = <IN>;
	close IN;

	@file = split(/\./, $file);
	$sample = $file[0];
	
	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($line =~ /^>/) {
			$count++;
		}
		$j++;
		$line=();
	}
	$j=0;
	
	print OUT "$sample\t$count\n";
	$count=0;
	$i++;
	$file=();
	@in=();
	@file=();
	$sample=();
}
$i=0;
close OUT;

}

####################

sub create_map_then_cat {

#var
my $i=0;
my $file;
my $sample;
my $filesize;
my $line;
my $readid;
my $j=0;

#array
my @output;
my @file;
my @in;

@output = qx(ls | grep fasta.trimmed);

open (OUT, ">>", "readid_sample.map") || die "Error can't open map file: $!\n";
open (OUT2, ">>", "cat.fasta") || die "Error can't open concatenated file: $!\n";

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	@file = split(/\./,$file);
	$sample = $file[0];

	$filesize = -s $file;

	if ($filesize > 0 ) {
		
		open (IN, "<", $file) || die "Error can't open fasta.trimmed file: $!\n";
		@in = <IN>;
		close IN;

		while ($in[$j]) {
			$line = $in[$j];
			chomp $line;

			if ($line =~ /^>/) {
				$line =~ /^>(\w{14,16})/;
				$readid = $1;
				print OUT "$readid\t$sample\n";
				print OUT2 "$line\n";
			}
			elsif (length($line) > 0) {
				print OUT2 "$line\n";
			}
			$j++;
			$line=();
			$readid=();
		}
		$j=0;
		@in=();
	}

	$i++;
	$file=();
	@file=();
	$sample=();
	$filesize=();
}
$i=0;
close OUT;
close OUT2;

}
