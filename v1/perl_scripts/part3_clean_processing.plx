#!/usr/bin/perl
#Dec.08,2010 by Terri Porter
#Merge post seqtrim steps into a single script to facilitate processing
#Need to have two files called *_seqtrim.fasta and *_seqtrim.qual in the cwd
#Usage $perl clean_processing.plx

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
uclust_clustering();
print "grab seed fastas to calculate stats\n";
grab_seed_fasta();
print "Calculating stats for seed fastas\n";
get_fasta_stats2();
print "filtering clusters for a minimum of 3 reads per cluster\n";
filter_clusters_by_read_frequency();
print "Calculating stats for filtered clusters\n";
get_fasta_stats3();
print "Mapping read number to cluster id\n";
create_cluster_read_map();
print "Starting blast searches, default local nt database, jobs spread to 8 cores\n";
blastall();
#submit jobs to megan for taxonomic summary and preliminary rarefaction

#######################

sub remove_newline_from_seqtrim_fasta {
	
	#var
	my $fasta;
	my $line;
	my $i=0;

	#array
	my @output;

	#@output = qx(ls | grep "seqtrim.fasta");
	#$fasta = $output[0];
	#chomp($fasta);

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

	#@output = qx(ls | grep "seqtrim.qual");
	#$qual = $output[0];
	#chomp($qual);

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

	@output = qx(ls | grep "seqtrim.fasta");
	$fasta = $output[0];
	chomp $fasta;

	open (FASTA,"<",$fasta) || die ("Error cannot open seqtrim fasta file to read: $!\n");
	#@fasta = <FASTA>;
	#close FASTA;

	@output2 = qx(ls | grep "seqtrim.qual");
	$qual = $output2[0];
	chomp $qual;

	open (QUAL,"<",$qual) || die ("Error cannot open seqtrim qual file to read: $!\n");
	#@qual = <QUAL>;
	#close QUAL;

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
		#$i++;
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
		#$j++;
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
	#@fasta = <FASTA>;
	#close FASTA;

	while (<FASTA>) {
		$line = $_;
		chomp($line);

		if ($flag==0) {
			if (/>/) {
				#$i++;
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
		#$i++;
	}

	push (@seq, $seq1); #don't forget to add last seq in file!

	while ($seq[$j]) {
		@split = split(//,$seq[$j]);
		$length = scalar(@split);
		push (@length, $length);
		@split = ();
		$j++;
	}
	
foreach (@length) {
	my $element = $_;
	print $element,"\n";
}

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

sub uclust_clustering {

	#array
	my @output;
	my @output2;

	@output = qx(usearch --cluster sorted.fna --uc results.uc --id 0.97 --nofastalign --rev --usersort --seedsout seed.fasta);
	#@output2 = qx(uclust --uc2fasta results.uc --input sorted.fna --output results.fasta);

}

############################

sub grab_seed_fasta {

	#var
	my $line;
	my $flag1 = 0;
	my $flag2 = 0;

	open (FASTA,"<","seed.fasta") || die ("Error cannot open results fasta file to read: $!\n");

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

	open (IN,"<","results_seed.fasta") || die ("Error cannot open results seed fasta file to read: $!\n");

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

	open (OUT, ">>", "results_seed.stats") || die ("Error cannot write to results seed stats file: $!\n");

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

sub blastall {

	#array
	my @output;

	@output = qx(blastall -p blastn -n T -d /net/bioinfo/1/blast/blastdb/nt -i sorted.fna.filtered -m 0 -o blastall.out -a 8);
	print "Blast search complete\n";
}

###########################

