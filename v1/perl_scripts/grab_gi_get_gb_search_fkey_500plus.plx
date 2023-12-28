#!/usr/bin/perl

#NGS_pipeline_part1.plx modified Feb.12,2013 and renamed as grab_gi_get_gb_search_fkey.plx

#Script to fetch sequences from GenBank using an Entrez query
#September 22, 2009 written by Terri Porter
#Usage $perl fetchseqs.plx

use warnings;
use strict;
use POSIX;

use Bio::DB::GenBank;     
use Bio::DB::Query::GenBank;
use Bio::SeqIO; 
use Bio::DB::EUtilities;
use Statistics::Lite qw(:all);

#perform entrez query, retrieve fasta
print "\nDoing entrez query...\n";
#var
my $query;
my $gb;
my $stream;
my $seqobj;
my $seq;
my $ID;
my $desc;

$query = Bio::DB::Query::GenBank->new         
	(-query   	=>'"Bovidae"[ORGN] AND rRNA[FKEY] AND "12S"',           
	-db      	=> 'nucleotide');

$gb = Bio::DB::GenBank->new
	(-format	=> 'Fasta');


$stream = $gb->get_Stream_by_query($query);     

open (OUT1,">>","entrez.fasta") || die ("Error cannot write to entrez.fasta: $!\n");

while( $seqobj =  $stream->next_seq() ) {       
 	$seq = $seqobj->seq();
	$ID = $seqobj->id();
	$desc = $seqobj->desc();
	print OUT1 ">".$ID.$desc."\n".$seq."\n";	
}     
close OUT1;

#parse fasta and get gi list and name
print "\nParsing fasta to get gi list...\n";
#var
my $i=0;
my $line;
my $gi;
#my $gb;
my $partial;
my $genus;
my $species;
my $gb_short;
my $gb_long;
my $name;

#array
my @in;
my @gi;
my @gb;
my @partial;
my @line;
my @gb_split;

#hash
my %gb_short_gi;
my %gb_short_gb_long;
my %gb_short_name;

open (IN,"<","entrez.fasta") || die ("Error cannot read from entrez.fasta: $!\n");
@in = <IN>;
close IN;

open (OUT2, ">>", "gi_gb.map") || die ("Error cannot write to gi_gb.map: $!\n");
#open (OUT2b,">>","gb.query") || die ("Error :$!\n");

open (OUT3,">>", "gi_genus_species.map") || die ("Error cannot write to gi_genus_species.map: $!\n");

while($in[$i]) {
	$line = $in[$i];
	chomp $line;
	if ($line =~ /^>/) {
		@line = split (/\|/,$line);
		$gi = $line[1];
		push(@gi,$gi);
		$gb = $line[3];
		push(@gb,$gb);
		$partial = $line[4];
		@partial = split(/ /,$partial);
		$genus = $partial[0];
		$species = $partial[1];
		$name = $genus." ".$species;

		#add hashes so I can track gb_short and gb_long and gi later
		@gb_split = split(/\./,$gb);
		$gb_short = $gb_split[0];
		$gb_short_gi{$gb_short} = $gi;
		$gb_short_gb_long{$gb_short} = $gb;
		print OUT2 "$gi\t$gb\n";
		print OUT3 "$gi\t$genus\t$species\n";
		$gb_short_name{$gb_short} = $name;
	}
	else {
		$i++;
		next;
	}
	$i++;
	$line=();
	@line=();
	$gi=();
	$gb=();
	$partial=();
	@partial=();
	$genus=();
	$species=();
	$name=();
	@gb_split=();
	$gb_short=();

}
close OUT2;
close OUT3;

#parse genbank records for features
print "\nGetting genbank records and searching for features...\n";
#var
my $factory;
my $file;
my $seqin;
my $seq2;
my $id;
my $feat_object;
my $sequence_fragment;
my $value;
my $SSU_seq;
my $ITS_seq;
my $LSU_seq;
my $spacer1;
my $spacer2;
my $scalar;
my $num_slices;
my $ceil;
my $k=1; #loop counter
my $batch=500;
my $temp;
my $max;
my $min;

#array
my @slicegi;

#get an error here if try to request too many at once so split into batches of 500
open (OUT4,">>","features.txt") || die ("Error cannot write to features.txt: $!\n");
print OUT4 "gb\tgi\tname\tsequence\n";

$scalar = scalar(@gi);
#print "scalar $scalar\n";
$num_slices = $scalar/500;
#print "num_slices $num_slices\n";
$ceil = ceil($num_slices);
#print "ceil $ceil\n";

#print "k $k\n";
while($k<=$ceil) {
	$temp = $k*$batch;
#	print "temp $temp\n";
	$max = $temp - 1;
#	print "max $max\n";
	$min = $max - 499;
#	print "min $min\n";
	@slicegi = @gi[$min..$max];
#	@slicegi = @gi[0..500];

#	print "@slicegi";

	grab_gb(\@slicegi,\$k); #pass reference to array
	
	$k++;
}
$k=0;

close OUT4;

##########
sub grab_gb {

my @slicegi = @{$_[0]}; #dereference
my $k = ${$_[1]}; #dereference

$factory = Bio::DB::EUtilities -> new (	-eutil	=>	'efetch',
					-db	=>	'nucleotide',
					-rettype =>	'gb',
					-id	=>	\@slicegi); #pass reference to array
$file = 'myseqs.gb'.$k;
$factory -> get_Response(-file => $file);
$seqin = Bio::SeqIO -> new (-file => $file,
			-format => 'genbank');

while ($seq2 = $seqin -> next_seq) {
	$id = $seq2 ->id; #looks like gb_short
	$gb_long = $gb_short_gb_long{$id}; #get associated gb_long
	$gi = $gb_short_gi{$id}; #get associated gi
	$name = $gb_short_name{$id}; #get associated name
	for $feat_object ($seq2->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "rRNA") {
			$sequence_fragment = $feat_object->spliced_seq->seq;
			if ($feat_object->has_tag('product')) {
				for $value ($feat_object -> get_tag_values('product')) {
					if ($value =~ /12S/) {
						$SSU_seq = $sequence_fragment; #keep var name as is, but it should probably be $12S_seq
					}
				}
			}
		}
	}
	print OUT4 "$gb_long\t$gi\t$name\t$SSU_seq\n";
	$SSU_seq="nil";
}

}
##########

#parse features.txt
print "\nGrabbing only 12S sequences that are at least 500bp...\n";
#var
my $j=0;
my $id2;
my $length1;
my $length2;
my $length_ITS;
my $length1b;
my $length1b_;
my $length2b;
my $length2b_;
my $min_length1;
my $min_length_ITS;
my $min_length2;
my $max_length1;
my $max_length_ITS;
my $max_length2;
my $mean_length1;
my $mean_length_ITS;
my $mean_length2;
my $min_length1b;
my $min_length_ITSb;
my $min_length2b;
my $max_length1b;
my $max_length_ITSb;
my $max_length2b;
my $mean_length1b;
my $mean_length_ITSb;
my $mean_length2b;

#array
my @file;
#my @line;
my @spacer1;
my @spacer2;
my @length1;
my @length2;
my @ITS_seq;
my @length_ITS;
my @spacer1b;
my @length1b;
my @spacer2b;
my @length2b;
my @ITS_seqb;
my @length_ITSb;

open (IN2,"<","features.txt") || die ("Error cannot read from features.txt: $!\n");
@file = <IN2>;
close IN2;

open (OUT5,">>","12S.fasta") || die ("Error cannot write to 12S.fasta :$!\n");

while ($file[$j]) {
	$line = $file[$j];
	chomp $line;
	@line = split(/\t/,$line);
	$gb_long = $line[0];
	$gi = $line[1];
	$name = $line[2];
	$SSU_seq = $line[3];
#	$spacer1 = $line[2];
#	$ITS_seq = $line[3];
#	$spacer2 = $line[4];
#	$LSU_seq = $line[5];
	if ($SSU_seq =~ /nil/) {
		$j++;
		next;
	}
#	elsif ($SSU_seq =~ /nil/) {
#		$j++;
#		next;
#		}
#	elsif ($ITS_seq =~/nil/) {
#		$j++;
#		next;
#	}
	else {
		@spacer1 = split(//,$SSU_seq); #should probably be @SSU_seq, but leave as @spacer1 for now
		$length1 = scalar(@spacer1);
		push(@length1,$length1);
#		@spacer2 = split(//,$spacer2);
#		$length2 = scalar(@spacer2);
#		push(@length2,$length2);
#		@ITS_seq = split(//,$ITS_seq);
#		$length_ITS = scalar(@ITS_seq);
#		push(@length_ITS,$length_ITS);
		if ($length1 >= 900) {
#			if ($length2 >= 100) {
#				@spacer1b = split(//,$spacer1);
#				$length1b = scalar(@spacer1b);
#				if ($length1b > 300) {
#					$length1b_ = $length1b."***";
#				}
#				else {
#					$length1b_ = $length1b;
#				}
#				push(@length1b,$length1b);
#				@spacer2b = split(//,$spacer2);
#				$length2b = scalar(@spacer2b);
#				if ($length2b > 300) {
#					$length2b_ = $length2b."***";
#				}
#				else {
#					$length2b_ = $length2b;
#				}
#				push(@length2b,$length2b);
#				@ITS_seqb = split(//,$ITS_seq);
#				$length_ITS = scalar(@ITS_seqb);
#				push(@length_ITSb,$length_ITS);
		print OUT5 ">$name $gb_long\n$SSU_seq\n";
#			}
		}
		else {
			$j++;
			next;
		}
	}
	$j++;
}
close OUT5;

#$min_length1 = min(@length1);
#$min_length_ITS = min(@length_ITS);
#$min_length2 = min(@length2);

#$max_length1 = max(@length1);
#$max_length_ITS = max(@length_ITS);
#$max_length2 = max(@length2);

#$mean_length1 = mean(@length1);
#$mean_length_ITS = mean(@length_ITS);
#$mean_length2 = mean(@length2);

#open (OUT6,">>","ITS.stats") || die ("Error cannot write to ITS.stats: $!\n");
#print OUT6 "ITS stats for all AFTOL sequences in Genbank with annotated spacer1, 5.8S, and spacer2 regions (before filtering for length)\n\n";
#print OUT6 "Feature\tMin\tMax\tMean\n";
#print OUT6 "ITS1:\t$min_length1\t$max_length1\t$mean_length1\n";
#print OUT6 "5.8S:\t$min_length_ITS\t$max_length_ITS\t$mean_length_ITS\n";
#print OUT6 "ITS2:\t$min_length2\t$max_length2\t$mean_length2\n";

#$min_length1b = min(@length1b);
#$min_length_ITSb = min(@length_ITSb);
#$min_length2b = min(@length2b);

#$max_length1b = max(@length1b);
#$max_length_ITSb = max(@length_ITSb);
#$max_length2b = max(@length2b);

#$mean_length1b = mean(@length1b);
#$mean_length_ITSb = mean(@length_ITSb);
#$mean_length2b = mean(@length2b);

#print OUT6 "\n\nafter filtering for length\n\n";
#print OUT6 "Feature\tMin\tMax\tMean\n";
#print OUT6 "ITS1:\t$min_length1b\t$max_length1b\t$mean_length1b\n";
#print OUT6 "5.8S:\t$min_length_ITSb\t$max_length_ITSb\t$mean_length_ITSb\n";
#print OUT6 "ITS2:\t$min_length2b\t$max_length2b\t$mean_length2b\n";

#close OUT6;

#resample fasta sequences for fragment analysis
#print "\nResampling fasta sequences for fragment analysis...\n";
#var
#my $k=0;
#my $l=0;
#my $m=0;
#my $fragment_size;
#my $outfile_name;
#my $header;
#my $sequence;
#my $base;
#my $cutoff;

#array
#my @fasta;
#my @fragment_sizes;
#my @sequence;
#my @sequence_fragment;
#my @sequence_reverse;
#my @sequence_fragment_reverse;

#open (IN3,"<","ITS.fasta") || die ("Error canot read from ITS.fasta: $!\n");
#@fasta = <IN3>;
#close IN3;

#foreach $line (@fasta) {
#	chomp $line;
#}

#@fragment_sizes = (50,100,150,200,250,300,350,400,450,500,550,600);

#create 5' fragments

#while ($fragment_sizes[$k]) {
#	$fragment_size = $fragment_sizes[$k];
	
#	$outfile_name = "5_prime_ITS_".$fragment_size.".fasta";
#	open (OUT7,">>",$outfile_name) || die ("Error cannot write to outfile: $!\n");
	
#	while ($fasta[$l]) {
#		$line = $fasta[$l];

#		if ($line =~ /^>/) {
#			$header = $line;
#			print OUT7 "$header\n";
#		}
#		else {
#			$sequence = $line;
#			@sequence = split(//,$sequence);
	
#			while ($sequence[$m]) {
#				$base = $sequence[$m];
#				$cutoff = $fragment_size-1;				
#				if ($m<=$cutoff) {
#					push(@sequence_fragment, $base);
#				}
#				$m++;
#			}
#			$m=0;
#			$sequence_fragment = join("",@sequence_fragment);
#			print OUT7 "$sequence_fragment\n";			
#			@sequence_fragment=();
#			@sequence=();
#		}
#		$l++;
#	}
#	close OUT7;
#	$l=0;
#	$k++;
#}
#$k=0;

#while ($fragment_sizes[$k]) {
#	$fragment_size = $fragment_sizes[$k];

#	$outfile_name = "3_prime_ITS_".$fragment_size.".fasta";
#	open (OUT8,">>",$outfile_name) || die ("Error cannot write to outfile2: $!\n");

#	while ($fasta[$l]) {
#		$line = $fasta[$l];

#		if ($line =~ /^>/) {
#			$header = $line;
#			print OUT8 "$header\n";
#		}
#		else {
#			$sequence = $line;
#			@sequence = split(//,$sequence);
#			@sequence_reverse = reverse(@sequence);

#			while ($sequence_reverse[$m]) {
#				$base = $sequence_reverse[$m];
#				$cutoff = $fragment_size-1;
#				if ($m<=$cutoff) {
#					push(@sequence_fragment_reverse, $base);
#				}
#				$m++;
#			}
#			$m=0;
#			@sequence_fragment = reverse(@sequence_fragment_reverse);
#			$sequence_fragment = join("",@sequence_fragment);
#			print OUT8 "$sequence_fragment\n";
#			@sequence_fragment_reverse=();
#			@sequence_fragment=();
#		}
#		$l++;
#	}
#	close OUT8;
#	$l=0;
#	$k++;
#}

