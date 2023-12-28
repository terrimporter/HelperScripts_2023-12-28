#!/usr/bin/perl
#Oct. 2, 2013 by Terri Porter
#Script to start with Shadi's illumina paired-end reads

use strict;
use warnings;

#global var
my $infile;
my $cluster = 1.0; #similarity cutoff for dereplication only
my $cores = 10; #for one cd-hit run at a time, or one usearch multithreaded job at a time
my $cluster_size_cutoff = 1; 

#global arrays
my @files;
my @simple_filenames;

#print "Get fasta files\n";
#get_fasta_files();

#print "Simplify filenames\n";
#simplify_filenames();

#print "Simplify readids, create readid_illuminareadid.map\n";
#simplify_readids(); #also create readid-illuminareadid map file 

#print "CD-HIT dereplication of each plate of reads\n";
#cluster(); #use cd-hit v.4.6i

#print "Creating clusterid map files\n";
#create_readid_cluster_map();

#print "Create fasta file formatted for usearch with cluster abundances\n";
#format_for_usearch();

print "Uchime denovo chimera filtering\n";
chimera_filter();

#print "Now, manually do taxonomic assignment using Insecta COI GenBank-genus trained NBC classifier.\n";

####################

sub get_fasta_files {

@simple_filenames = qx(ls | grep fasta);

#print "@files\n";

}

####################

sub simplify_filenames {

my $i=0;
my $file;
my @file;
my $plate;
my $newfilename;
my $output;

while ($files[$i]) {
	$file = $files[$i];
	chomp $file;

	@file = split(/_/,$file);
	$plate = $file[0];

	$newfilename = $plate.".fasta";

	push (@simple_filenames, $newfilename);

	$output = qx(mv $file $newfilename);

	$i++;

}
$i=0;

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

	$newfilename = $file.".simplified";

	open (OUT, ">>", $newfilename) || die "Error cannot open new simplified file: $!\n";

	open (MAP, ">>", "readid_illuminareadid.map") || die "Error cannot open readid - illumina readid map file: $!\n";

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($line =~ /^>/) {
			$line =~ /^>(\S+)/;
			$illumina_readid = $1;
			$simple_readid = $readid_counter;
			
			print MAP "$simple_readid\t$illumina_readid\n";
			
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
	@in=();
	$newfilename=();
	$file=();
	close OUT;
	close MAP;
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

my @output;
my $i=0;
my $file;
my $outfile;
my $output;

@output = qx(ls | grep simplified);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	$outfile = $file.".derep";

	$output = qx(/home/terri/cd-hit-v4.6-2012-04-25/cd-hit-est -i $file -o $outfile -c $cluster -n 8 -M 53000 -T $cores -r 0);
	$i++;
}
$i=0;

}

####################

sub create_readid_cluster_map {

my @output;
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
my $k=0;
my $file;

@output = qx(ls | grep clstr);

while ($output[$k]) {
	$file = $output[$k];
	chomp $file;

	open (IN, "<", $file) || die "Error cannot open cd hit clstr file:$!\n";
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

	$k++;
	$file=();
	@in=();
}
$k=0;

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
my $regex;
my $k=0;
my $file;
my $outfile;
my @output;

open (MAP, "<", "clusterreadid_readabundance.map") || die "Error cannot open clusterreadid-readabundance.map: $!\n";
@map=<MAP>;
close MAP;

$regex = "'derep"."\$'";
@output = qx(ls | grep $regex);

while ($output[$k]) {
	$file = $output[$k];
	chomp $file;

	open (FASTA, "<", $file) || die "Error cannot open $file: $!\n";
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
	
	$outfile = $file.".abund";
	open (OUT, ">>", $outfile) || die "Error cannot open $outfile: $!\n";

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
	close OUT;
	$k++;
	$file=();
	@fasta=();
}
$k=0;

}

####################

sub chimera_filter {

my $output;
my @output;
my $i=0;
my $file;
my $outfile;
my $nonch_outfile;
my $ch_outfile;
my $prefix;
my $outfile2;

@output = qx(ls | grep '.abund\$');
#print "output:@output\n";

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;
	
	if ($file =~ /Plate\d+_/) {
		$file =~ /(Plate\d+)_/;
		$prefix = $1;
#		print "prefix:$prefix\n";

		$outfile = $prefix.".sorted";
		print "outfile1:$outfile\n";

		$output = system("usearch6.0.307_i86linux32 -sortbysize $file -output $outfile");

		$nonch_outfile = $prefix.".nonch.fasta";
		$ch_outfile = $prefix.".ch.fasta";

		$output = system("usearch6.0.307_i86linux32 -uchime_denovo $outfile -nonchimeras $nonch_outfile -chimeras $ch_outfile");

		$outfile2 = $prefix."_derep_minsize".$cluster_size_cutoff.".fasta";
		$output = system("usearch6.0.307_i86linux32 -sortbysize $nonch_outfile -minsize $cluster_size_cutoff -output $outfile2");

	}
	
	$i++;
	$file=();
	$prefix=();
	$outfile=();
	$output=();
	$nonch_outfile=();
	$ch_outfile=();
	$outfile2=();

}
$i=0;

}

####################
