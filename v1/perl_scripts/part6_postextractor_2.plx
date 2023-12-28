#!/usr/bin/perl
#Feb. 1, 2011 edited to grab the whole ITS from the start of ITS1 for MEGAN classification
#Dec. 6, 2011 by Terri Porter
#Script to process files after running FungalITSextractor
#Ensure ITS1.fasta, outdata.csv, and seqtrim.qual.newlineremoved are available in date/
#usage perl part6_postextractor.plx seqtrim.fasta.newlineremoved outdata.csv seqtrim.qual.newlineremoved

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $ref;
my $outfile;

#declare array
my @FASTA;
my @outdata;
my @qual;
my @ITS_filtered;
my @ITS_reformatted;
my @qual_trimmed;

#declare hash

open (FASTA, "<", $ARGV[0]) || die "Error cannot open seqtrim.fasta.newlineremoved: $!\n";
@FASTA = <FASTA>;
close FASTA;

open (OUTDATA, "<", $ARGV[1]) || die "Error cannot open outdata.csv: $!\n";
@outdata = <OUTDATA>;
close OUTDATA;

open (QUAL, "<", $ARGV[2]) || die "Error cannot open seqtrim.qual.newlineremoved: $!\n";
@qual = <QUAL>;
close QUAL;

$outfile = "seqtrim.fasta.stats";
print "running get_stats\n";
get_stats(\@FASTA,$outfile);

print "running length_filter\n";
length_filter(\@FASTA);

open (ITS_FILTERED, "<", "ITS.filtered") || die "Error cannot open ITS.filtered: $!\n";
@ITS_filtered = <ITS_FILTERED>;
close ITS_FILTERED;

$outfile = "ITS.filtered.stats";
print "running get_stats2\n";
get_stats(\@ITS_filtered);

print "running trim_qual_by_outdata\n";
trim_qual_by_outdata(\@ITS_filtered,\@outdata,\@qual);

print "running reformat_cluster_header3\n";
reformat_cluster_header3(\@ITS_filtered);

open (ITS_REFORMATTED, "<", "ITS.filtered.reformatted") || die "Error cannot read from ITS.filtered.reformatted: $!\n";
@ITS_reformatted = <ITS_REFORMATTED>;
close ITS_REFORMATTED;

open (QUAL_TRIMMED,"<","seqtrim.qual.newlineremoved.trimmed") || die "Error cannot read from seqtrim.qual.newlineremoved.trimmed: $!\n";
@qual_trimmed = <QUAL_TRIMMED>;
close QUAL_TRIMMED;

print "running consolidate\n";
consolidate(\@ITS_reformatted,\@qual_trimmed);

####################

sub consolidate {

#declare var
my $i=0;
my $line;
my $header;
my $id;
my $j;
my $seq;
my $qual;
my $mean;
my $num_elements;
my $line1;
my $line2;
my $x;
my $num;
my $check=();
my $length;
my $array_ref1 = $_[0];
my $array_ref2 = $_[1];

#declare array
my @fasta = @$array_ref1;
my @qual = @$array_ref2;
my @sorted_length;
my @sorted_qual;
my @phred;
my @unsorted;
my @sorted;
my @x;

#declare hash
my %fasta;
my %qual;
my %mean;
my %num;
my %map;
my %map_sorted;

#put fasta file into hash indexed by id
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w{14})/;
		$id = $1;
		$j=$i+1;
		$seq = $fasta[$j];
		chomp $seq;
		$fasta{$id} = $seq;
	}
	$i++;
}
$i=0;

#put qual file into hash indexed by id
while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w{14})/;
		$id = $1;
		$j=$i+1;
		$qual = $qual[$j];
		chomp $qual;
		$qual{$id} = $qual;

		@phred = split(/\s+/,$qual);
		$num_elements = scalar (@phred);
		$mean = mean (@phred);
		$mean{$id} = $mean;
		$num{$id} = $num_elements;
	}
	$i++;
}
$i=0;

#create array that includes id,length,avephred
while (($id,$num) = each (%num) ) {
	$mean = $mean{$id};
	$line = $id.",".$num.",".$mean;
	push(@unsorted,$line);
}

#Schwartzian Transform to sort multiple columns
@sorted = map { $_->[0]}
sort { $b->[2] <=> $a->[2] || $b->[3] <=> $a->[3] } 
map { [$_, split(/,/)] }
@unsorted;

open (FASTA,">>","file.fasta") || die "Error cannot write to fasta.file $!\n";
open (QUAL,">>","file.qual") || die "Error cannot write to fasta.qual $!\n";

foreach $x (@sorted) {
		@x = split(/,/,$x);
		$id = $x[0];
		if ($fasta{$id}) {
			$length = $num{$id};
			if ($length >= 50) { ##### only keep ITS2 if longer than 50 bp #####
				$seq = $fasta{$id};
				print FASTA ">$id\n$seq\n";
				$qual = $qual{$id};
				print QUAL ">$id\n$qual\n";
			}
		}
		$length=();
}
print "-----> printed file.fasta and file.qual\n";

}

####################

sub reformat_cluster_header3 {

#declare var
my $line;
my $clusterID;
my $symbol;
my $rest;
my $id;
my $array_ref = $_[0];
my $i=0;

#declare array
my @line;
my @input = @$array_ref;

open (OUT,">>","ITS.filtered.reformatted") || die "Error cannot write to ITS.filtered.reformatted: $!\n";

while ($input[$i]) {
	$line = $input[$i];
	chomp $line;
	if ($line =~ /^>/) {
		if ($line =~ />\w{14}.+/) {
			#@line = split(/\|/,$line);
			#$clusterID = $line[2];
			#$symbol = $line[3];
			#$rest = $line[4];
			$line =~ />(\w{14}).+/;
			$id = $1;
			print OUT ">$id\n";
		}
	}
	else {
		print OUT "$line\n";
	}
	$i++;
}
$i=0;
close OUT;
print "-----> printed ITS.filtered.reformatted\n";

}

####################

sub trim_qual_by_outdata {

#declare var
my $i=0;
my $line;
my $header;
my $id;
my $ITS1_trim;
my $start="nil";
my $end="nil";
my $id_qual;
my $j;
my $seq;
my $tot_length;
my $start_adj;
my $end_adj;
my $k=0;
my $base;
my $array_ref1 = $_[0];
my $array_ref2 = $_[1];
my $array_ref3 = $_[2];
my $trimmed_seq;
my $id_fasta;

#declare array
my @ids;
my @line;
my @seq;
my @fasta = @$array_ref1;
my @outdata = @$array_ref2;
my @qual = @$array_ref3;
my @trimmed_seq;

#declare hash
my %id_start;
my %id_end;

#grab list of ids from ITS.filtered
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>\w+/) {
		$header = $line;
		$line =~ /^>(\w+)/;
		$id = $1;
#		$id = substr $header, 1, 14;
		push (@ids,$id);
	}
	$i++;
}
$i=0;

#parse trim data from the whole outdata file
while ($outdata[$i]) {
	$line = $outdata[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$id = $line[0];
	$id =~ /\(.+\)\s+(\w+)/;
	$id = $1;
	$id = substr $id, 0, 14;	
	$ITS1_trim = $line[5];
	if ($ITS1_trim !~ /-----/) {
		$ITS1_trim =~ /ITS1:\s+(\d+)-(end|\d+)/;
		$start = $1;
		$end = "end"; # instead of $end = $2;
		#print "start: $start\t end: $end\n";#test
	}
	else {
		$start = "nil";
		$end = "nil";
	}
	$id_start{$id} = $start;
	$id_end{$id} = $end;
	#test
	#print "id: $id\t trim_start: $start\t trim_end: $end\n";
	$i++;
}
$i=0;

open (OUT,">>","seqtrim.qual.newlineremoved.trimmed") || die "Error cannot write to seqtrim.qual.newlineremoved.trimmed: $!\n";

#for each fasta id, grab trim info from hashes, edit seqtrim.qual
while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w+)\s+/;
		$id_qual = $1;
		#print "id_qual: $id_qual\n";#test
		$j = $i+1;
		$seq = $qual[$j];
		chomp $seq;
		#print "seq: $seq\n";#test
		@seq = split(/\s+/,$seq);
		$tot_length = scalar(@seq);
		#print "tot_length: $tot_length\n";#test
		#print "seq: @seq\n";#test

		$start = $id_start{$id_qual};
		#print "start: $start\n";#test
		if ($start eq "nil") {
			$i++;
			next;
		}
		else {
			$start_adj = $start-1;
			#print "start_adj: $start_adj\n";#test
		}

		$end = $id_end{$id_qual};
		if ($end eq "end") {
			$end = $tot_length;
			$end_adj = $end-1;
		}
		elsif ($end eq "nil") {
			$i++;
			next;
		}
		else {	
			$end_adj = $end-1;
		}

		while ($seq[$k]) {
			$base = $seq[$k];
			#print "base: $base\t k: $k\t start_adj: $start_adj\t end_adj: $end_adj\n";
			#print "$base ";
			if ($k >= $start_adj) {
				#print "$k >= $start_adj\n";
				if ($k <= $end_adj) {
					#print "$k<= $end_adj\n";
					push(@trimmed_seq,$base);
					#print "@trimmed_seq\n";
				}
				else {
					$k++;
					next;
				}
			}
			else {
				$k++;
				next;
			}
			$k++;
		}
		#print "@trimmed_seq\n";#test
		#print "\n";
		$k=0;
		$trimmed_seq = join(" ",@trimmed_seq);
		#print "trimmed_seq: $trimmed_seq\n";#test
		print OUT ">$id_qual\n$trimmed_seq\n";
	}
	$i++;
	@trimmed_seq=();
	$id_qual=();
	$seq=();
	$start="nil";
	$end="nil";
	$base=();
	$trimmed_seq=();
}
$i=0;
close OUT;
print "-----> printed seqtrim.qual.newlineremoved.trimmed\n";

open (OUT2,">>","seqtrim.fasta.newlineremoved.trimmed") || die "Error cannot write to seqtrim.fasta.newlineremoved.trimmed: $!\n";

#for each fasta id, grab trim info from hashes, edit seqtrim.fasta
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w+)\s+/;
		$id_fasta = $1;
		#print "id_qual: $id_qual\n";#test
		$j = $i+1;
		$seq = $fasta[$j];
		chomp $seq;
		#print "seq: $seq\n";#test
		@seq = split(//,$seq);
		$tot_length = scalar(@seq);
		#print "tot_length: $tot_length\n";#test
		#print "seq: @seq\n";#test

		$start = $id_start{$id_fasta};
		#print "start: $start\n";#test
		if ($start eq "nil") {
			$i++;
			next;
		}
		else {
			$start_adj = $start-1;
			#print "start_adj: $start_adj\n";#test
		}

		$end = $id_end{$id_fasta};
		if ($end eq "end") {
			$end = $tot_length;
			$end_adj = $end-1;
		}
		elsif ($end eq "nil") {
			$i++;
			next;
		}
		else {	
			$end_adj = $end-1;
		}

		while ($seq[$k]) {
			$base = $seq[$k];
			#print "base: $base\t k: $k\t start_adj: $start_adj\t end_adj: $end_adj\n";
			#print "$base ";
			if ($k >= $start_adj) {
				#print "$k >= $start_adj\n";
				if ($k <= $end_adj) {
					#print "$k<= $end_adj\n";
					push(@trimmed_seq,$base);
					#print "@trimmed_seq\n";
				}
				else {
					$k++;
					next;
				}
			}
			else {
				$k++;
				next;
			}
			$k++;
		}
		#print "@trimmed_seq\n";#test
		#print "\n";
		$k=0;
		$trimmed_seq = join(" ",@trimmed_seq);
		#print "trimmed_seq: $trimmed_seq\n";#test
		print OUT2 ">$id_fasta\n$trimmed_seq\n";
	}
	$i++;
	@trimmed_seq=();
	$id_fasta=();
	$seq=();
	$start="nil";
	$end="nil";
	$base=();
	$trimmed_seq=();
}
$i=0;
close OUT2;
print "-----> printed seqtrim.fasta.newlineremoved.trimmed\n";


}

####################

sub length_filter {

#declare var
my $i=0;
my $line;
my $header;
my $j;
my $next_line;
my $length;
my $sequence;
my $min_length=50; ##### set minimum sequence length here #####
my $array_ref = $_[0];

#declare array
my @next_line;
my @input = @$array_ref;

open (OUT,">>","ITS.filtered") || die ("Error creating outfile: $!\n");

while ($input[$i]) {
	$line = $input[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$j=$i+1;
		$next_line = $input[$j];
		chomp $next_line;
		@next_line = split(//,$next_line);
		$length = scalar(@next_line);

		if ($length >= $min_length) {
			$sequence = join('',@next_line);
			print OUT "$header\n$sequence\n";
			$sequence=();
		}
		$i+=2;
		@next_line=();
	}
	else {
		$i++;
	}
}
close OUT;
print "-----> printed ITS.filtered\n";

}

####################

sub get_stats {

#declare variables
my $line;
my $flag=0;
my $seq1;
my $seq2;
my $seq;
my $i=0;
my $length;
my $min;
my $max;
my $mean;
my $mode;
my $array_ref = $_[0];

#declare array
my @seq;
my @split;
my @length;
my @input = @$array_ref;

while($input[$i]){
	$line = $input[$i];
	chomp $line;

	if ($flag==0){

		if ($line =~ />/){
			$i++;
			next;
		}
		else {
			$seq1 = $line;
			$flag=1;
		}	
	}	
	elsif ($flag==1){
		
		if ($line =~ />/){
			$flag=0;
			push (@seq, $seq1);
		}
		else {
			$seq2 = $line;
			$seq1 = $seq1.$seq2; ##process non-header files this way to accomodate incorrectly formatted fasta files
		}
	}
	$i++;
}
push (@seq, $seq1);#don't forget to add last seq in file!
$i=0;

while ($seq[$i]){
	@split = split(//, $seq[$i]);
	$length = scalar(@split);
	push (@length, $length);
	@split =();#empty array
	$i++;
}
$i=0;

$min = min (@length);
$max = max (@length);
$mean = mean (@length);
$mode = mode (@length);
my $num = scalar(@seq);

print "-----> printing stats\n";

open (OUT,">>",$outfile) || die "Error cannot print fasta stats outfile: $!\n";

print OUT "NumSeqs\tMin\tMax\tMean\tMode\n";
print OUT $num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";

close OUT;
}

####################
