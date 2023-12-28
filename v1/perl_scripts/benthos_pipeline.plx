#!/usr/bin/perl
#Oct. 2, 2013 by Terri Porter
#Script to start with Shadi's illumina paired-end reads, already paired with PANDA and qc'd with printseq
#create illumina readid -> simple readid map file
#create simple readid -> sample map file
#merge all
#cd hit cluster 98%
#create simple readid -> cluster readid map file
#remove singletons
#usearch sorty by size, uchime_denovo chimera filtering
#NBC GenBank-genus training set, classify usign 70% cutoff (for 200bp frags assigned to genus)
#MEGAN LCA parse
#export cluster readid-> taxonid
#grab lineage
#map back cluster readids to lineage
#for each genus, grab one rep seq, multiply by number of sites
#track abundance of OTUs for map file
#track sample for map file
#align using MAFFT
#remove ambiguously aligned characters using trimAl
#tree using RAxML with GTR + gamma
#UNIFRAC needs newick tree file with matching map file with samples and abundances

#/1/scratch/terri/Biomonitoring_benthos_100113/extracted_fasta/
#USAGE perl benthos_pipeline.plx

use strict;
use warnings;

#global var
my $infile;
my $cluster = 0.98;
my $cores = 15;

#global arrays
my @files;
my @simple_filenames;

print "Get fasta files\n";
get_fasta_files();

print "Simplify filenames\n";
simplify_filenames();

print "Simplify readids, create 2 mapping files\n";
simplify_readids(); #also create readid-illuminareadid map file AND readid-sample map file

print "Merge fastas\n";
merge_fastas();

print "Cluster fastas\n";
$infile = "cat.fasta";
cluster(\$infile); #use cd-hit v.4.6i

print "Creating clusterid map files\n";
create_readid_cluster_map();

print "Create fasta file formatted for usearch with cluster abundances\n";
format_for_usearch();

print "Uchime denovo chimera filtering\n";
chimera_filter();

print "Now, manually do taxonomic assignment using Insecta COI GenBank-genus trained NBC classifier.\n";

####################

sub get_fasta_files {

@files = qx(ls | grep fasta);

#print "@files\n";

}

####################

sub simplify_filenames {

my $i=0;
my $file;
my $site;
my $locus;
my $subsite;
my $newfilename;
my $output;

while ($files[$i]) {
	$file = $files[$i];
	chomp $file;

	$file =~ /(PAD|WC)(-|_)(AD|BE)(-|_)(\d{1,2}(A|B|C))/;
	$site = $1;
	$locus = $3;
	$subsite = $5;
	$newfilename = $site."_".$locus."_".$subsite.".fasta";

	push (@simple_filenames, $newfilename);

	$output = qx(mv $file $newfilename);

#	print "$site\t$locus\t$subsite\n";

	$i++;

}
$i=0;
#print "@simple_filenames\n";

}

####################

sub simplify_readids {

my $i=0;
my $file;
my $site;
my $locus;
my $subsite;
my $j=0;
my $newfilename;
my $line;
my $illumina_readid;
my $readid_counter = 1;
my $simple_readid;

my @in;

while ($simple_filenames[$i]) {
	$file = $simple_filenames[$i];
	chomp $file;

	open (IN, "<", $file) || die "Error cannot open simple file name: $!\n";
	@in = <IN>;
	close IN;

	$file =~ /(PAD|WC)_(AD|BE)_(\d{1,2}(A|B|C))/;
	$site = $1;
	$locus = $2;
	$subsite = $3;

	$newfilename = $site."_".$locus."_".$subsite.".simplified";

	open (OUT, ">>", $newfilename) || die "Error cannot open new simplified file: $!\n";

	open (MAP, ">>", "readid_illuminareadid.map") || die "Error cannot open readid - illumina readid map file: $!\n";

	open (MAP2, ">>", "readid_sample.map") || die "Erorr cannot open readid - sample map file: $!\n";

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($line =~ /^>/) {
			$line =~ /^>(\S+)/;
			$illumina_readid = $1;
			$simple_readid = $readid_counter;
			
			print MAP "$simple_readid\t$illumina_readid\n";
			print MAP2 "$simple_readid\t$site$subsite\n";
			
			$readid_counter++;
			$line = ">$simple_readid\n";
			print OUT "$line";
		}
		else {
			print OUT "$line\n";
		}
		$j++;
	}
	$j=0;
#	$readid_counter=1;
	@in=();
	$newfilename=();
	$file=();
	close OUT;
	close MAP;
	close MAP2;
	$i++;
}
$i=0;

}

####################

sub merge_fastas {

my $i=0;
my $file;
my @in;
my $line;
my $j=0;

@files=();

@files = qx(ls | grep simplified);

open (OUT, ">>", "cat.fasta") || die "Error cannot open cat.fasta: $!\n";

while ($files[$i]) {
	$file = $files[$i];
	chomp $file;

	open (IN, "<", $file) || die "Error cannot open simplified infile: $!\n";
	@in = <IN>;
	close IN;

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if (length($line) > 0) {
			print OUT $line."\n";
		}
		$j++;
	}
	$j=0;
	$i++;
	$file=();
	@in=();
}
$i=0;
close OUT;

}

####################

sub cluster {

$infile = ${$_[0]};
my $output;

$output = qx(/home/terri/cd-hit-v4.6-2012-04-25/cd-hit-est -i $infile -o cluster.fasta -c $cluster -n 8 -M 53000 -T $cores -r 0);

}

####################

sub create_readid_cluster_map {

my @in;
my $i=0;
my $line;
my %set=(); #reads that belong to OTU, indexed by readid
my @line;
my $readid;
my @readids;
my $readid_line;
my $original;
my $new;
my $secondpartline;
my $cluster_readid;
my $scalar;
my $count;
my $j;
my $flag=0;

open (IN, "<", "cluster.fasta.clstr") || die "Error cannot open cd hit clstr file:$!\n";
@in = <IN>;
close IN;

open (OUT2, ">>", "clusterreadid_readabundance.map") || die "Error cannot open clusterreadid - readabundance: $!\n";

open (OUT3, ">>", "clusterreadid_readid.map") || die "Error cannot open clusterreadid - readid:$!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0 && $line =~ /^>/) {
		$flag=1;
		$i++;
		next;
	}

	if ($flag==1 && $line !~ /^>/) {
		if ($line =~ /\*/) {
			@line = split(/\t/,$line);
			$secondpartline = $line[1];
			$secondpartline =~ />(\d+)\.+/;
			$cluster_readid = $1;
			$set{"clusterreadid"} = $cluster_readid;	
		}
		elsif ($line !~ /\*/) {
			@line = split(/\t/,$line);
			$secondpartline = $line[1];
			$secondpartline =~ />(\d+)\.+/;
			$readid = $1;
			if (exists $set{"readid"}) {
				$original = $set{"readid"};
				$new = $original."|".$readid;
				$set{"readid"} = $new;
			}
			else {
				$set{"readid"} = $readid;
			}
		}
		$i++;
		next;
	}
	elsif ($flag==1 && $line =~ /^>/) {
	
		#process previous full set		
		if (exists $set{"clusterreadid"}) { #don't forget to parse last set
			$cluster_readid = $set{"clusterreadid"};
			print OUT3 "$cluster_readid\t$cluster_readid\n";
		}

		if (exists $set{"readid"}) {
			$readid_line = $set{"readid"};
			@readids = split(/\|/,$readid_line);
			$scalar = scalar(@readids);
			foreach $readid (@readids) {
				print OUT3 "$cluster_readid\t$readid\n";
			}
			$count = $scalar+1;
			print OUT2 "$cluster_readid\t$count\n";
		}
		else {
			print OUT2 "$cluster_readid\t1\n";
		}
			
		%set=();
		@readids=();
		$scalar=();
		$count=();
		
		$flag=0;
		$i--;
	}
	$i++;
	$line=();
}
$i=0;

#process the last full set
if (exists $set{"clusterreadid"}) { #don't forget to parse last set
	$cluster_readid = $set{"clusterreadid"};
	print OUT3 "$cluster_readid\t$cluster_readid\n";
}

if (exists $set{"readid"}) {
	$readid_line = $set{"readid"};
	@readids = split(/\|/,$readid_line);
	$scalar = scalar(@readids);
	foreach $readid (@readids) {
		print OUT3 "$cluster_readid\t$readid\n";
	}
	$count = $scalar+1;
	print OUT2 "$cluster_readid\t$count\n";
}

close OUT2;
close OUT3;

}

####################

sub format_for_usearch {

my @map;
my @fasta;
my $i=0;
my $line;
my @line;
my $cluster_readid;
my $abundance;
my %abundance; #indexed by cluster_readid
my $j;
my $nextline;
my $newline;

open (MAP, "<", "clusterreadid_readabundance.map") || die "Error cannot open clusterreadid-readabundance.map: $!\n";
@map=<MAP>;
close MAP;

open (FASTA, "<", "cluster.fasta") || die "Error cannot open cluster.fasta: $!\n";
@fasta = <FASTA>;
close FASTA;

#hash map file
while ($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$cluster_readid = $line[0];
	$abundance = $line[1];
	$abundance{$cluster_readid} = $abundance;
	$i++;
}
$i=0;

open (OUT, ">>", "cluster.fasta.abund") || die "Error cannot open cluster.fasta.abund: $!\n";

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\d+)/;
		$cluster_readid = $1;

		if (exists $abundance{$cluster_readid}) {
			$abundance = $abundance{$cluster_readid};
			
			if ($abundance == 0) { #change to ==1 to remove singletons
				$i+=2;
				next;
			}
			else {
				$newline = $line.";size=".$abundance.";";
				print OUT "$newline\n";
				$j = $i+1;
				$nextline = $fasta[$j];
				chomp $nextline;
				print OUT "$nextline\n";
			}
		}
		else {
			print "Erorr: cannot find abundance for $cluster_readid\n";
		}
	}
	$i++;
	$line=();
	$cluster_readid=();
	$abundance=();
	$nextline=();
}
$i=0;

}

####################

sub chimera_filter {

my $output;

$output = system('usearch6.0.307_i86linux32 -sortbysize cluster.fasta.abund -output cluster.fasta.abund.sorted');

$output = system('usearch6.0.307_i86linux32 -uchime_denovo cluster.fasta.abund.sorted -nonchimeras cluster.nonch.fasta -chimeras cluster.ch.fasta');

$output = system('usearch6.0.307_i86linux32 -sortbysize cluster.nonch.fasta -minsize 2 -output cluster_98otus_minsize2.fasta');

}

####################

