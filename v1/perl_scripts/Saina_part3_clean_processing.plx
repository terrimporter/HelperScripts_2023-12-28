#!/usr/bin/perl
#Dec.08,2010 by Terri Porter
#Merge post seqtrim steps into a single script to facilitate processing
#Need to have two files called *_seqtrim.fasta and *_seqtrim.qual in the cwd
#Usage $perl clean_processing.plx
#modify to work with Saina_part1_sort.plx

use strict;
use warnings;

#global var

print "starting to remove empty fastas from seqtrim\n";
remove_empty_fastas_from_seqtrim_output();
print "starting to remove improper newline endings\n";
remove_newline_from_seqtrim_fasta();
remove_newline_from_seqtrim_qual();
print "calculating seqtrim fasta stats\n";
get_fasta_stats();
print "sorting by quality for uclust input\n";
sort_by_qual();
print "starting uclust dereplication/clustering step\n";
usearch_clustering();
#print "grab seed fastas to calculate stats\n";
#grab_seed_fasta();
print "Calculating stats for seed fastas\n";
get_fasta_stats2();
print "filtering clusters for a minimum of 3 reads per cluster\n";
filter_clusters_by_read_frequency();
print "Calculating stats for filtered clusters\n";
get_fasta_stats3();
print "Mapping read number to cluster id\n";
create_cluster_read_map();
print "Starting megablast searches, Saina's database\n";
megablast();
print "Starting to parse megablast tabular output\n";
parse_megablast();
print "Compiling hit_read_frequencies\n";
hit_read_frequency();
print "Reformatting hit read frequency table\n";
full_hit_read_frequency();

#######################

sub remove_newline_from_seqtrim_fasta {
	
	#var
	my $fasta;
	my $line;
	my $i=0;

	#array
	my @output;

	open (FASTA,"<","seqtrim.fasta.noempty") || die ("Error cannot open seqtrim fasta noempty file to read: $!\n");

	open (TEMP,">>","seqtrim.fasta.newlineremoved.temp") || die ("Error cannot open seqtrim fasta newlineremoved temp file to write: $!\n");

	while (<FASTA>) {
		$line = $_;
		chomp($line);
		if ($line =~ /^>/) {
			print TEMP "\n",$line,"\n";
		}
		else {
			print TEMP $line;
		}
	}
	close FASTA;
	close TEMP;

	open (IN,"<","seqtrim.fasta.newlineremoved.temp") || die ("Error cannot open seqtrim fasta newlineremoved temp file to read: $!\n");
	
	open (OUT,">>","seqtrim.fasta.newlineremoved") || die ("Error cannot open seqtrim fasta newlineremoved file to write: $!\n");

	while (<IN>) {
		$line = $_;
		chomp($line);

		if ($i==0) {
			$i=1;
			next;
		}
		else {
			print OUT $line."\n";
		}
	}
	close IN;
	close OUT;

	unlink("seqtrim.fasta.newlineremoved.temp");
}

#####################

sub remove_newline_from_seqtrim_qual {

	#var
	my $qual;
	my $line;
	my $flag=0;
	my $i=0;

	#array
	my @output;

	open (QUAL,"<","seqtrim.qual.noempty") || die ("Error cannot open seqtrim qual file to read: $!\n");

	open (TEMP,">>", "seqtrim.qual.newlineremoved.temp") || die ("Error cannot open seqtrim qual newlineremoved temp file to write: $!\n");

	while (<QUAL>) {
		$line = $_;
		chomp($line);
		if ($line =~ /^>/) {
			print TEMP "\n",$line,"\n";
			$flag=0;
		}
		elsif ($flag==0) {
			print TEMP $line;
			$flag=1;
		}
		elsif ($flag>0) {
			print TEMP " $line";
		}
	}
	close QUAL;
	close TEMP;

	open (IN,"<","seqtrim.qual.newlineremoved.temp") || die ("Error cannot open seqtrim qual newlineremoved temp file to read: $!\n");

	open (OUT,">>","seqtrim.qual.newlineremoved") || die ("Error cannot open seqtrim qual newlineremoved file to write: $!\n");

	while (<IN>) {
		$line = $_;
		chomp($line);

		if ($i==0) {
			$i++;
			next;
		}
		else {
			print OUT $line."\n";
		}
	}
	close IN;
	close OUT;

	unlink("seqtrim.qual.newlineremoved.temp");

}

#####################

sub remove_empty_fastas_from_seqtrim_output {
	
	#var
	my $fasta;
	my $qual;
	my $line;
	my $i=0;
	my $j=0;

	#array
	my @output;
	my @output2;
	my @fasta;
	my @qual;
	my @line;

	@output = qx(ls | grep "seqtrim");
	$fasta = $output[0];
	$qual = $output[1];
	chomp $fasta;
	chomp $qual;

	open (FASTA,"<",$fasta) || die ("Error cannot open seqtrim fasta file to read: $!\n");

	open (QUAL,"<",$qual) || die ("Error cannot open seqtrim qual file to read: $!\n");

	while (<FASTA>) {
		$line = $_;
		chomp($line);
		if ($line =~ /^>/) {
			push (@line, $line);
		}
		else {
			if ($line =~ /\w+/) {
				push (@line, $line);
			}
			else {
				pop (@line);
			}
		}
	}
	close FASTA;

	open (OUT1,">>","seqtrim.fasta.noempty") || die ("Error cannot open seqtrim fasta noempty file: $!\n");

	foreach (@line) {
		$line = $_;
		print OUT1 $line."\n";
	}
	close OUT1;
	
	@line=();

	while (<QUAL>) {
		$line = $_;
		chomp($line);
		if ($line =~ /^>/) {
			push (@line, $line);
		}
		else {
			if ($line =~ /\w+/) {
				push (@line, $line);
			}
			else {
				pop (@line);
			}
		}
	}
	close QUAL;

	open (OUT2,">>","seqtrim.qual.noempty") || die ("Error cannot open seqtrim qual noempty file: $!\n");

	foreach (@line) {
		$line = $_;
		print OUT2 $line."\n";
	}
	close OUT2;
}

################################

sub get_fasta_stats {

	use Statistics::Lite qw(:all);

	#var
	my $line;
	my $flag=0;
	my $seq1;
	my $seq2;
	my $i=0;
	my $j=0;
	my $length;
	my $element;
	my $min;
	my $max;
	my $mean;
	my $mode;
	my $num;

	#array
	my @fasta;
	my @seq;
	my @split;
	my @length;

	open (FASTA,"<","seqtrim.fasta.newlineremoved") || die ("Error cannot open seqtrim fasta newlineremoved file to read: $!\n");

	while (<FASTA>) {
		$line = $_;
		chomp($line);

		if ($flag==0) {
			if (/>/) {
				next;
			}
			else {
				$seq1 = $line;
				$flag=1;
			}
		}
		elsif ($flag==1) {
			if (/^>/) {
				$flag=0;
				push (@seq,$seq1);
			}
			else {
				$seq2 = $line;
				$seq1 = $seq1.$seq2;
			}
		}
	}

	push (@seq, $seq1); #don't forget to add last seq in file!

	while ($seq[$j]) {
		@split = split(//,$seq[$j]);
		$length = scalar(@split);
		push (@length, $length);
		@split = ();
		$j++;
	}
	
#foreach (@length) {
#	my $element = $_;
	#print $element,"\n";
#}

	open (OUT,">>","seqtrim.fasta.newlineremoved.stats") || die ("Error cannot open seqtrim fasta newlineremoved stats file to write: $!\n");

	$min = min (@length);
	$max = max (@length);
	$mean = mean (@length);
	$mode = mode (@length);
	$num = scalar (@seq);

	print OUT "NumSeqs\tMin\tMax\tMean\tMode\n";
	print OUT $num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";
}

#############################

sub sort_by_qual {

	#var
	my $i=0;
	my $line;
	my $j=0;
	my $phred_score;
	my $k=0;
	my $value;
	my $l;
	my $index;
	my $header;
	my $phred_seq;
	my $m=0;
	my $fna_header;
	my $fna_seq;

	#array
	my @fna;
	my @qual;
	my @header;
	my @phred_seq;
	my @phred_scores;
	my @value;
	my @fna_header;
	my @fna_seq;

	#hash
	my %hash;

	open (FASTA,"<","seqtrim.fasta.newlineremoved") || die ("Error cannot open seqtrim fasta newlineremoved file: $!\n");
	@fna = <FASTA>;
	close FASTA;
	#test
	my $scalar = scalar(@fna);
	my $numlines = $scalar/2;
	print "num lines in fasta: $numlines\n"; #test

	open (QUAL,"<","seqtrim.qual.newlineremoved") || die ("Error cannot open seqtrim qual newlineremoved file: $!\n");
	@qual = <QUAL>;
	close QUAL;

	while ($qual[$i]) {
		$line = $qual[$i];
		chomp $line;
		if ($line =~ /^>/) {
			push (@header,$line);
		}
		else {
			push (@phred_seq, $line);
		}
		$i++;
	}

	while ($phred_seq[$j]){
		$line = $phred_seq[$j];
		@phred_scores = split (/ /,$line);
		foreach $phred_score (@phred_scores) {
			if ($phred_score < 20 ) {
				$k++;
			}
		}
		$hash{$j} = $k;
		@phred_scores=();
		$k=0;
		$j++;
	}

	foreach $value (sort {$hash{$a} <=> $hash{$b}} keys %hash ) {
	       push (@value, $value);
     	}

	open (OUT1,">>","sorted.qual") || die ("Error cannot open sorted qual file: $!\n");

	foreach $l (@value) {
		$index = $l;
		$header = $header[$index];
		$phred_seq = $phred_seq[$index];
		print OUT1 "$header\n$phred_seq\n";
	}
	close OUT1;

	while ($fna[$m]) {
		$line = $fna[$m];
		chomp $line;
		if ($line =~ /^>/) {
			push (@fna_header, $line);
		}
		else {
			push(@fna_seq, $line);
		}
		$m++;
	}

	open (OUT2,">>","sorted.fna") || die ("Error cannot open sorted fna file: $!\n");

	foreach my $l (@value) {
		$index = $l;
		$fna_header = $fna_header[$index];
		$fna_seq = $fna_seq[$index];
		print OUT2 "$fna_header\n$fna_seq\n";
	}
	close OUT2;
}

################################

sub usearch_clustering {

	#array
	my @output;
	my @output2;

	@output = qx(usearch --cluster sorted.fna --uc results.uc --id 1.0 --nofastalign --rev --usersort --seedsout seed.fasta);
	#@output2 = qx(usearch --uc2fastax results.uc --input sorted.fna --output results.fasta);

}

############################

sub grab_seed_fasta {

	#var
	my $line;
	my $flag1 = 0;
	my $flag2 = 0;

	open (FASTA,"<","results.fasta") || die ("Error cannot open results fasta file to read: $!\n");

	open (OUT,">>","results_seed.fasta") || die ("Error cannot open results seed fasta file to write: $!\n");

	while (<FASTA>) {
		$line = $_;
		chomp($line);
		if ($line =~ /\*/) {
			print OUT "$line\n";
			$flag1 = 1;
			$flag2 = 0;
		}
		elsif ($flag1 == 1) {
			print OUT "$line\n";
			$flag1=0;
			$flag2=1;
		}
		elsif ($flag2==1) {
			next;
		}
	}
	close FASTA;
	close OUT;
}

############################

sub get_fasta_stats2 {

	#var
	my $line;
	my $flag=0;
	my $seq1;
	my $seq2;
	my $i=0;
	my $length;
	my $element;
	my $min;
	my $max;
	my $mean;
	my $mode;
	my $num;

	#array
	my @seq;
	my @split;
	my @length;

	open (IN,"<","seed.fasta") || die ("Error cannot open results seed fasta file to read: $!\n");

	while (<IN>) {
		$line = $_;
		chomp($line);

		if ($flag==0) {
			if (/>/) {
				next;
			}
			else {
				$seq1 = $line;
				$flag=1;
			}
		}
		elsif ($flag==1) {
			if (/>/) {
				$flag=0;
				push (@seq,$seq1);
			}
			else {
				$seq2 = $line;
				$seq1 = $seq1.$seq2;
			}
		}
	}
	push (@seq, $seq1);
	close IN;

	while ($seq[$i]) {
		@split = split (//,$seq[$i]);
		$length = scalar (@split);
		push (@length, $length);
		@split = ();
		$i++;
	}

	foreach (@length) {
		$element = $_;
		print $element, "\n";
	}

	$min = min (@length);
	$max = max (@length);
	$mean = mean (@length);
	$mode = mode (@length);
	$num = scalar (@seq);

	open (OUT, ">>", "seed.stats") || die ("Error cannot write to results seed stats file: $!\n");

	print OUT "NumSeqs\tMin\tMax\tMean\tMode\n";
	print OUT $num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";
}

#############################

sub filter_clusters_by_read_frequency {

	#var
	my $line;
	my $cluster_size;
	my $id_line;
	my $id;
	my $j=0;
	my $i=0;

	#array
	my @line;
	my @ids;

	open (IN1,"<","results.uc") || die ("Error cannot read from results uc file: $!\n");

	while (<IN1>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^C/) {
			@line = split (/\t/,$line);
			$cluster_size = $line[2];
			$id_line = $line[8];
			if ($cluster_size >= 3) { ###EDIT READ FREQUENCY CUTOFF HERE###
				$id_line =~ /^(\w{14,16})\s+/;
				$id = $1;
				push (@ids, $id);
			}
		}
	}
	close IN1;

	open (IN2,"<","sorted.fna") || die ("Error cannot read from sorted fna file: $!\n");

	open (OUT, ">>", "sorted.fna.filtered") || die ("Error cannot write to sorted fna filtered file: $!\n");

	while (<IN2>) {
		$line = $_;
		chomp $line;
		if ($ids[$j]){
			if ($i==0) {
				$id = $ids[$j];
				if ($line =~ /$id/) {
					print OUT $line."\n";
					$i=1;
					next;
				}
				else {
					$i=0;
					next;
				}
			}
			elsif ($i==1) {
				print OUT $line."\n";
				$i=0;
			}
			$j++;
		}
	}
	close IN2;
	close OUT;

}

##############################

sub get_fasta_stats3 {

	use Statistics::Lite qw(:all);

	#var
	my $line;
	my $flag=0;
	my $seq1;
	my $seq2;
	my $i=0;
	my $length;
	my $element;
	my $num;
	my $min;
	my $max;
	my $mean;
	my $mode;

	#array
	my @seq;
	my @split;
	my @length;

	open (IN,"<","sorted.fna.filtered") || die ("Error cannot read from sorted fna filtered file: $!\n");

	open (OUT,">>","sorted.fna.filtered.stats") || die ("Error cannot write to sorted fna filtered stats file: $!\n");

	while (<IN>) {
		$line = $_;
		chomp($line);

		if ($flag==0) {

			if (/>/) {
				next;
			}
			else {
				$seq1 = $line;
				$flag=1;
			}
		}
		elsif ($flag==1) {
			
			if (/>/) {
				$flag=0;
				push (@seq, $seq1);
			}
			else {
				$seq2 = $line;
				$seq1 = $seq1.$seq2;
			}
		}
	}
	close IN;
	push (@seq, $seq1);

	while ($seq[$i]) {
		@split = split (//,$seq[$i]);
		$length = scalar(@split);
		push (@length, $length);
		@split = ();
		$i++;
	}

	foreach (@length){
		$element = $_;
		print OUT $element,"\n";
	}

	$min = min(@length);
	$max = max(@length);
	$mean = mean(@length);
	$mode = mode(@length);
	$num = scalar(@seq);

	print OUT "NumSeqs\tMin\tMax\tMean\tMode\n";
	print OUT $num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";
	close OUT;
	
}	
#############################

sub create_cluster_read_map {
	
	#var
	my $line;
	my $id_line;
	my $cluster_size;
	my $id;

	#array
	my @line;

	open (IN,"<","results.uc") || die ("Error cannot read from results uc file: $!\n");

	open (OUT,">>","results.uc.clusters") || die ("Error cannot read from results uc clusters file: $!\n");

	print OUT "ID\tCluster_size\n";

	while (<IN>) {
		$line = $_;
		chomp($line);
		if ($line =~ /^C/) {
			@line = split (/\t/,$line);
			$id_line = $line[8];
			#print "id_line test: $id_line\n";#test
			my @id_line = split (/ /,$id_line);
			$cluster_size = $line[2];
			$id = $id_line[0];
			#print "test id: $id\n";#test
			print OUT "$id\t$cluster_size\n";
		}
		else {
			next;
		}
	}
	close IN;
	close OUT;

}

##############################

sub megablast {

	#array
	my @output;

	@output = qx(blastn -task megablast -db /home/terri/blast-2.2.24_32bit/db/Saina -query sorted.fna.filtered -out megablast.out -outfmt 6);
	print "Megablast complete.\n";
}

###########################

sub parse_megablast {

	#var
	my $i=0;
	my $line;
	my $identity;
	my $queryID_current;
	my $queryID_previous = "nil";
	my $flag=0;
	my $key;
	my $value;
	my $j=0;
	my $bitscore_entry;
	my $values_index;
	my $line_index;
	my $best_line;
	my $best_bitscore;
	my $x;
	my $size;

	#declare array
	my @file;
	my @line;
	my @temp;
	my @values;
	my @keys;
	my @entry;

	#declare hashes
	my %bitscore;
	my %line;

	open (IN, "<", "megablast.out") || die ("Error cannot read megablast.out: $!\n");
	@file = <IN>;
	close IN;

	open (OUT, ">>", "megablast.filtered") || die ("Error cannot write to megablast.filtered: $!\n");

	#filter for %identities >= 99%
	while ($file[$i]) {
		$line = $file[$i];
		chomp $line;
		@line = split (/\t/,$line);
		$identity = $line[2];
		if ($identity >= 99) { #####modify percent identity cutoff here#####
			push (@temp,$line);
		}
		$i++;
	}

	#for each query, filter by bitscore
	$i=0;
	while ($temp[$i]) {
		$line = $temp[$i];
		$line{$i} = $line;
		@entry = split(/\t/,$line);
		$queryID_current = $entry[0];
		$bitscore_entry = $entry[11];
		if ($queryID_current ne $queryID_previous) {
			if ($flag==0) {
				$bitscore{$i} = $bitscore_entry;
				$queryID_previous = $queryID_current;
				$flag=1;
			}
			elsif ($flag==1) {
				foreach $key (sort{$bitscore{$b}<=>$bitscore{$a}} keys %bitscore) {
					$value = $bitscore{$key};
					push (@values,$value);
					push (@keys,$key);
				}
				$best_bitscore = $values[0];
				if ($best_bitscore >= 100) {##### modify bitscore cutoff here #####
					$values_index = $keys[0];
					$x = $values_index-$i;
					$line_index = $i + $x;
					$best_line = $line{$line_index};
					print OUT "$best_line\n";
				}
				@values=();
				@keys=();
				%bitscore=();
				$bitscore{$i}=$bitscore_entry;
				$queryID_previous=$queryID_current;
			}
		}
		elsif ($queryID_current eq $queryID_previous) {
			$bitscore{$i}=$bitscore_entry;
			$queryID_previous = $queryID_current;
		}
		$i++;
	}

	#dont' forget to parse last set of entries!
	foreach $key (sort{$bitscore{$b}<=>$bitscore{$a}} keys %bitscore) {
		$value = $bitscore{$key};
		push(@values,$value);
		push(@keys,$key);
	}
	$best_bitscore = $values[0];
	if ($best_bitscore >= 100) {##### modify bitscore cutoff here too #####
		$values_index = $keys[0];
		$x = $values_index-$i;
		$line_index = $i+$x;
		$best_line = $line{$line_index};
		print OUT "$best_line\n";
	}
	close OUT;
}

####################

sub hit_read_frequency {

	#var
	my $line;
	my $i=0;
	my $current_id;
	my $id_line;
	my $id;
	my $ref;
	my $cluster_size;
	my $current_ref;
	my $k;
	my $j=0;
	my $total_frequency;

	#array
	my @line;
	my @ids;
	my @refs;
	my @refs_nonredundant;
	my @cluster_sizes;

	#hash
	my %refs;

	#get filtered list of ids and refs

	open (IN1, "<", "megablast.filtered") || die ("Error cannot read megablast.filtered: $!\n");

	while (<IN1>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^G/) {
			@line = split (/\t/,$line);
			$id_line = $line[0];
			$id_line =~ /^(\S{14})/;
			$id = $1;
			push (@ids,$id);
			$ref = $line[1];
			push (@refs, $ref);
		}
		else {
			next;
		}
	}
	close IN1;
	print "array ids: @ids\n";#test
	print "array refs: @refs\n";#test

	#get list of refs and cluster_sizes

	open (TMP,">>","temp.txt") || die ("Error cannot open temp file: $!\n");

	while ($ids[$i]) {
		$current_id = $ids[$i];
		open (IN2,"<","results.uc.clusters") || die ("Error cannot open results.uc.clusters: $!\n");

		while (<IN2>) {
			$line = $_;
			chomp $line;
			if ($line =~ /$current_id/) {
				@line = split (/\t/,$line);
				$cluster_size = $line[1];
				my $refvar = $refs[$i];
				print TMP "$refvar\t$cluster_size\n";
			}
			else {
				next;
			}
		}
		close IN2;
		$i++;
	}
	close TMP;
	
	#dereplicate @refs
	
	%refs = map {$_,1} @refs;
	@refs_nonredundant = keys %refs;

	#print hit read frequency table
		
	open (OUT,">>","hit_read_freq.txt") || die ("Error can't write to hit_read_freq.txt: $!\n");
	print OUT "ReferenceID\tReadFrequency\n";

	while ($refs_nonredundant[$j]) {
		$current_ref = $refs_nonredundant[$j];
		open (IN3,"<","temp.txt") || die ("Error can't write to temp.txt: $!\n");
			
		while (<IN3>) {
			$line = $_;
			chomp $line;
			if ($line =~ /$current_ref/) {
				@line = split (/\t/,$line);
				$cluster_size = $line[1];
				push (@cluster_sizes,$cluster_size);
			}
			else {
				next;
			}
		}
		$j++;
		$total_frequency = 0;
		foreach $k (@cluster_sizes) {
			$total_frequency += $k;
		}
		print OUT "$current_ref\t$total_frequency\n";
		@cluster_sizes=();
	}
	close IN3;
	close OUT;
	unlink("temp.txt");
}

####################

sub full_hit_read_frequency {

	#var
	my $reference_fasta;
	my $line;
	my $ref_ID;
	my $i=0;
	my $current_ref;
	my $read_freq;
	my $total;

	#array
	my @ref_IDS;
	my @line;
	my @read_freq;

	print "Please enter name of reference.fasta:\n";
	$reference_fasta = <STDIN>;
	chomp $reference_fasta;

	open (IN1,"<", $reference_fasta) || die ("Error cannot read reference.fasta: $!\n");

	while (<IN1>) {
		$line = $_;
		chomp $line;
		if ($line =~ />\S+/) {
			$line =~ />(\S+)/;
			$ref_ID = $1;
			my @ids=split(/\|/,$ref_ID);
			my $part_ID = $ids[0];
			push (@ref_IDS, $part_ID);
		}
		else {
			next;
		}
	}
	print "array ref IDS: @ref_IDS\n";#test

	close IN1;

	open (OUT,">>", "full_hit_read_freq.table") || die ("Error cannot write to full hit read freq table: $!\n");
	print OUT "ReferenceID\tReadFrequency\n";

	while ($ref_IDS[$i]) {
		$current_ref = $ref_IDS[$i];
		print "var current_ref: $current_ref\t";#test
		open (IN2, "<", "hit_read_freq.txt") || die ("Error: $!\n");

		while (<IN2>) {
			$line = $_;
			chomp $line;
			if ($line =~ /$current_ref/) {
				@line = split (/\t/,$line);
				$read_freq = $line[1];
				print "var read freq: $read_freq\n";#test
				push(@read_freq,$read_freq);
			}
			else {
				$read_freq=0;
				push(@read_freq,$read_freq);
				next;
			}
		}
		close IN2;
		foreach (@read_freq){
			$total += $_;
		}
		print OUT "$current_ref\t$total\n";
		$total=0;
		$i++;
		@read_freq=();
	}

}
