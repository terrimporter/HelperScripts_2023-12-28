#!/usr/bin/perl
# FASTBLAST_xml.plx
#Terri Porter, May 2004
#updated March 2006 because GenBank changed the format of their text blast reports
#updated August 2009 to retrieve xml formatted blast output to be retrieved later
#updated January 2010 updated to put all outfiles into a separate directory, prints a hit table to summarize results, prints a list of gi numbers that can be used with get_genbank_record.plx and parse_genbank_records_using_gi.plx to get lineage info
#updated March 2010 to grab sequences from top blast hits

#SCRIPT TO BLASTX A FASTA FORMATTED TEXT FILE AGAINST GENBANK

use strict;
use warnings;
use lib "/projects/cipres/bioperl-1.4";
use Bio::Tools::Run::RemoteBlast;
use Bio::Seq;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Bio::SearchIO;
use Bio::SeqIO;
use IO::String;

#declare some variables
my $count = 0;
#counts number of sequences to BLAST
#my $i = 1; 
#keeps track of accession numbers in array
my $v = 1; 
#verbose
my $accession;
#declare global variables
my $path;
my $path2;

#declare arrays
my @accession;
my @new_array;

#requires arguments to specify the sequence infile
my $infile= $ARGV[0]; 
#holds the filename of the input file
unless (@ARGV) { 
	#holds arguments
	print "You need to enter a command line argument like this:
	perl <script> <inputfile>\n";
	exit;
}

#define new SeqIO object
my $seqio = Bio::SeqIO -> new (
	-file => $infile,
	-format => 'fasta')
or die "Could not create Bio::SeqIO\n";

#add these paramaters to standard blast factory object
$Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{'ALIGNMENTS'} = '1';
$Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{'HITLIST_SIZE'} = '1';
$Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{'FORMAT_TYPE'} = 'XML';
$Bio::Tools::Run::RemoteBlast::HEADER{'SERVICE'} = 'megablast';

#create remote blast factory object and initialize blast parameters
my $blast_factory = Bio::Tools::Run::RemoteBlast -> new (
	'-prog' => 'blastn',
	'-data' => 'nr',
	'-expect' => '10.0',
	'-readmethod' => 'xml');

#make a separate directory to hold the blast reports
mkdir ("fastblast", 0755) || print $!;

#loop over all sequences input and count how many sequences are found 
while (my $seq = $seqio -> next_seq) {
	$count++;
	#loop to blast each sequence, in turn, against the database
	my $job = $blast_factory -> submit_blast ($seq);
	print STDERR "Blasting sequence number $count\n";
	#loop to load rids returned for the blast job submitted
	while (my @rids = $blast_factory -> each_rid) {
		#loop over rids to check server for a result
		foreach my $rid (@rids) {
			my $blast_results = $blast_factory -> retrieve_blast($rid);
			#print "$blast_results\n";
			if ($blast_results == 0) { 
				#still waiting to complete search
				print STDERR "."; #watch dots while waiting
				sleep 5; #pause between checking for results
			}
			elsif ($blast_results == (-1)) {
				#error returned, remove RID from stack
				print STDERR "retrieve_blast returns -1\n";
				$blast_factory -> remove_rid($rid);
			}
			#use Bio::SearchIO to return a Bio::SearchIO result object
			else {
				print STDERR "Receiving blast results...\n";
				$blast_results -> verbose (0);
				my $result = $blast_results -> next_result();
				my $filename = $result -> query_name()."\.OUT";
				$path = "/home/biodept/tp45/fastblast";
				chdir($path) || die "$!";
				$blast_factory -> save_output ($filename);
				my $path2 = "/home/biodept/tp45";
				chdir($path2) || die "$!";
				$blast_factory -> remove_rid($rid);
			}
		}
	}
}

chdir($path) || die "$!";

#subroutine to merge blast reports into single file called merge.txt
merge();
sub merge {
	my $dir = "$path";
	opendir DH, $dir;
	my @files = readdir (DH);

	#open merged file in append mode
	open (OUT, ">>merge.txt") ||die ("Cannot open new merged file");

	foreach my $files (@files) {
		#open the files one by one
		open (FH, $files) || die ("Cannot open file in array");
		#read the contents
		@new_array = <FH>;
		#store in the merged file
		print OUT @new_array;
		close (FH);
	}
	close (OUT);
}

if (@new_array) {
	#print "\nBlast reports have been successfully merged into a single file.\n\n";
}

#subroutine to parse merge.txt and produce three outfiles: 
#1) hit.table
#2) ids.infile which contains a list of gi numbers
#list of ids can be used with get_genbankrecord.plx and parse_genbankrecords_using_gi.plx to get lineage info
#3) hit.fasta
parse_top_xml_hits();
sub parse_top_xml_hits {
	my $i=0;
	open(IN, "<merge.txt");
        open(OUT, ">>hit.table");
	open(OUT2, ">>ids.infile");
	open (OUT3, ">>hit.fasta");
	print OUT "ContigID\tgi\tHitDef\tHspAcc\tHspBitScore\tHspEval\tHspQueryFrom\tHspQueryTo\tHspHitFrom\tHspHitTo\tHspIdentity\tHspAlnLen\n";
	while (<IN>) {
		my ($line) = $_;
		chomp ($line);
		if ($i==0){		#ContigID
			if ($line =~ /<Iteration_query-def>(.+)<\/Iteration_query-def>/) {
				my($querydef) = $1;
				if ($querydef =~ /(.+)/) {
					my($defline) = $1;
					#my($length) = $2;
					#my($numreads) = $3;
					print OUT "$defline\t";
					$i=1;
					next;
				}
			}
		}
		if ($i==1) {		#gi
			if ($line =~ /<Hit_id>(.+)<\/Hit_id>/) {
				my($HitId) = $1; 
				if ($HitId =~ /gi\|(\d+)\|gb\|(\S{10})\|/) {
					my($gi) = $1;
					print OUT "$gi\t";
					print OUT2 "$gi\n";
					$i=2;
					next;
				}
			}
		}
		if ($i==2) {			#HitDef
			if ($line =~ /<Hit_def>(.*)<\/Hit_def>/) {
				my($hitdef) = $1;		
				if ($hitdef =~ /(.*)/) {
					my($def) = $1;
					#my($species) = $2;
					print OUT "$def\t";
					$i=3;
					next;
				}
			}
		}
		if ($i==3) {			#HitAcc
			if ($line =~ /<Hit_accession>(.*)<\/Hit_accession>/) {
				my($accession) = $1;
				print OUT "$accession\t";
				print OUT3 ">$accession\n";
				$i=4;
				next;
			}
		}
		if ($i==4) {			#HspBitScore
			if ($line =~ /<Hsp_bit-score>(.*)<\/Hsp_bit-score>/) {
				my($bitscore) = $1;
				print OUT "$bitscore\t";
				$i=5;
				next;
			}
		}
		if ($i==5) {			#HspEval
			if ($line =~ /<Hsp_evalue>(.*)<\/Hsp_evalue>/) {
				my($evalue) = $1;
				print OUT "$evalue\t";
				$i=6;
				next;
			}
		}
		if ($i==6) {			#HspQueryFrom
			if ($line =~ /<Hsp_query-from>(.*)<\/Hsp_query-from>/) {
				my($HspQueryFrom) = $1;
				print OUT "$HspQueryFrom\t";
				$i=7;
				next;
			}
		}
		if ($i==7) {			#HspQueryTo
			if ($line =~ /<Hsp_query-to>(.*)<\/Hsp_query-to>/) {
				my($HspQueryTo) = $1;
				print OUT "$HspQueryTo\t";
				$i=8;
				next;
			}
		}
		if ($i==8) {			#HspHitFrom
			if ($line =~ /<Hsp_hit-from>(.*)<\/Hsp_hit-from>/) {
				my($HspHitFrom) = $1;
				print OUT "$HspHitFrom\t";
				$i=9;
				next;
			}
		}
		if ($i==9) {			#HspHitTo
			if ($line =~ /<Hsp_hit-to>(.*)<\/Hsp_hit-to>/) {
				my($HspHitTo) = $1;
				print OUT "$HspHitTo\t";
				$i=10;
				next;
			}
		}
		if ($i==10) {			#HspIdentity
			if ($line =~ /<Hsp_identity>(.*)<\/Hsp_identity>/) {
				my($HspIdentity) = $1;
				print OUT "$HspIdentity\t";
				$i=11;
				next;
			}
		}	
		if ($i==11) {			#HspAlnLen
			if ($line =~ /<Hsp_align-len>(.*)<\/Hsp_align-len>/) {
				my($HspAlignLen) = $1;
				print OUT "$HspAlignLen\n";
				$i=12;
				next;
			}
		}
		if ($i==12) {			#HspHitSeq
			if ($line =~ /<Hsp_hseq>(.*)<\/Hsp_hseq>/) {
				my($HspHitSeq) = $1;
				print OUT3 "$HspHitSeq\n";
				$i=0;
			}
		}
	}
        print "Outfiles can be found in the folder fastblast/\n\n";
}

