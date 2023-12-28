#!/usr/bin/perl
#Dec.21,2010 by Terri Porter
#Script to reformat uclust results.fasta to a .list file to be read by mothur
#need to put seqtrim.fna and seqtrim.qual cleaned sequences into working directory, make sure these files have empties and newlines removed
#usage $perl make_list_file.plx results.fasta

use strict;
use warnings;

#declare global array
my @output;

#declare global var
my $basename;

print "Getting list of infiles names.\n";
get_names_of_infiles();
print "List of infiles complete.\n";

print "Make .groups file for mothur.\n";
make_groups_file();
print ".groups file written.\n";

print "Merging files.\n";
merge_files();
print "Merge complete.\n";

print "Start sorting reads by quality.\n";
sort_reads_by_qual();
print "Read quality sorting complete.\n";

print "Starting clustering.\n";
uclust();
print "Clustering complete.\n";

print "Making .list file.\n";
make_list_file();
print ".list file written.\n";
print "\nReady to manually use .group and .list files with mothur.\n";
####################

sub get_names_of_infiles {

	#declare var
	my $seqtrim;
	
	@output = qx(ls | grep ".fna.seqtrim");
	foreach $seqtrim (@output) {
		chomp $seqtrim;
	}
}

####################

sub make_groups_file {

	#declare var
	my $i=0;
	my $seqtrim;
	my $line;
	my $id;
	my $group;
	my $groups_file;
	
	#declare array
	my @seqtrim;
	my @line;
	
	print "Please enter base name for .groups and .list files (ex. 16S):\n";
	$basename = <STDIN>;
	chomp $basename;
	$groups_file = $basename.".groups";
	
	open (OUT,">",$groups_file) || die ("Error cannot write to group file: $!\n");
	
	while ($output[$i]) {
		$seqtrim = $output[$i];
		@seqtrim = split (/\./,$seqtrim);
		$group = $seqtrim[0];
		
		open (IN,"<",$seqtrim) || die ("Error cannot read seqtrim file:$!\n");
		
		while (<IN>) {
			$line = $_;
			chomp $line;
			if (/>/) {
				@line = split(/ /,$line);
				$id = $line[0];
				$id =~ s/>//;
				print OUT "$id\t$group\n";
			}
		}
		close IN;
		$i++;
	}
	close OUT;
}

####################

sub merge_files {

	#declare var
	my $pwd;
	my $dir;
	my $x;
	my $filename;
	my $pathtofile;
	my $dir_name;
	my $oldfilepath;
	my $newfilepath;

	#declare array
	my @output2;
	my @output3;
	my @output4;
	my @files;
	my @fna;
	my @qual;
	my @new_array;
	my @new_array2;
	
	use File::Copy;
	
	$dir_name = "to_merge";
	mkdir $dir_name;
		
	@output2 = qx(pwd);
	$pwd = $output2[0];
	chomp $pwd;
	
	print "pwd: $pwd\n";#test
	
	$dir = $pwd."/to_merge/";

	print "dir: $dir\n";#test
	
	@output3 = qx(ls | grep "seqtrim");
	
	print "output3: @output3\n";#test
	
	foreach $x (@output3) {
		chomp $x;
		$oldfilepath = $pwd."/".$x;
		$newfilepath = $dir."/".$x;
		move ($oldfilepath, $newfilepath);
	}

	opendir DH, $dir;
	@files = readdir(DH);
	
	print "files: @files\n";#test
	
	foreach $x (@files) {
		if ($x=~ /\.fna/) {
			push (@fna, $x);
		}
		elsif ($x =~ /\.qual/) {
			push (@qual, $x);
		}
	}
	
	print "fna: @fna\n";#test
	print "qual: @qual\n";#test
	
	open (OUT1, ">>","merge.fna.txt") || die ("Cannot write to merge.fna.txt: $!\n");
	
	foreach $x (@fna) {
		$filename = $x;
		$pathtofile = $dir."/".$filename;
		open (FH, $pathtofile) || die ("Cannot read file in to_merge: $!\n");
		@new_array = <FH>;
		print OUT1 @new_array;
		close FH;
	}
	close OUT1;

	open (OUT2, ">>", "merge.qual.txt") || die ("Cannot write to merge.qual.txt: $!\n");

	foreach $x (@qual) {
		$filename = $x;
		$pathtofile = $dir."/".$filename;
		open (FH2, $pathtofile) || die ("Cannot read file in to_merge: $!\n");
		@new_array2 = <FH2>;
		print OUT2 @new_array2;
		close FH2;
	}
	close OUT2;

}

###################

sub sort_reads_by_qual {

	#declare var
	my $i=0;
	my $line;
	my $j=0;
	my $test;
	my $test2;
	my $test3;
	my $test4;
	my $test5;
	my $phred_score;
	my $k=0;
	my $l;
	my $index;
	my $header;
	my $phred_seq;
	my $m=0;
	my $fna_header;
	my $fna_seq;
	my $key;
	my $value;
	
	#declare array
	my @fna;
	my @qual;
	my @header;
	my @phred_seq;
	my @phred_scores;
	my @phred_count;
	my @value;
	my @fna_header;
	my @fna_seq;
	
	#declare hash
	my %hash;
	
	open (FNA,"<","merge.fna.txt") || die ("Error cannot read from merge.fna.txt: $!\n");
	@fna = <FNA>;
	close FNA;

	open (QUAL,"<","merge.qual.txt") || die ("Error cannot read from merge.qual.txt: $!\n");
	@qual = <QUAL>;
	close QUAL;

	while ($qual[$i]) {
		$line = $qual[$i];
		chomp $line;
		if ($line =~ /^>/) {
			push (@header, $line);
		}
		else {
			push (@phred_seq,$line);
		}
		$i++;
	}

	$test = scalar(@header);
	$test2 = scalar(@phred_seq);
	print "header: $test\nphred_seq: $test2\n";

	while ($phred_seq[$j]) {
		$line = $phred_seq[$j];
		@phred_scores = split (/ /, $line);
		foreach $phred_score (@phred_scores) {
			if ($phred_score < 20) {
				$k++;
			}
		}

		$hash{$j} = $k;
		@phred_scores =();
		$k=0;
		$j++;
	}

	foreach $value (sort {$hash{$a} <=> $hash{$b}} keys %hash) {
		push (@value, $value);
	}
	
	$test3 = scalar(@value);
	print "value: $test3\n";

	open (QUAL2,">>", "sorted.qual") || die ("Error cannot write to sorted.qual: $!\n");
	
	foreach my $l (@value) {
		$index = $l;
		$header = $header[$index];
		$phred_seq = $phred_seq[$index];
		print QUAL2 "$header\n$phred_seq\n";
	}
	close QUAL2;
	
	while ($fna[$m]) {
		$line = $fna[$m];
		chomp $line;
		if ($line =~ /^>/) {
			push (@fna_header, $line);
		}
		else {
			push (@fna_seq, $line);
		}
		$m++;
	}

	$test4 = scalar(@fna_header);
	$test5 = scalar(@fna_seq);
	print "fna_header: $test4\nfna_seq: $test5\n";

	open (FNA2, ">>", "sorted.fna") || die ("Error cannot write to sorted.fna: $!\n");

	foreach $l (@value) {
		$index = $l;
		$fna_header = $fna_header[$index];
		$fna_seq = $fna_seq[$index];
		print FNA2 "$fna_header\n$fna_seq\n";
	}
	close FNA2;

}

####################

sub uclust {

	#declare array
	my @output2;
	my @output3;

	@output2 = qx(usearch --cluster sorted.fna --uc results.uc --id 0.97 --nofastalign --rev --usersort);
	@output3 = qx(usearch --uc2clstr results.uc --input sorted.fna --output results.clstr);
	
}

####################
sub make_list_file {
	
#declare var
my $results_file;
my $label;
my $num_otus;
my $line;
my $id_part;
my $id;
my $group;
my $outfile1;
my $outfile2;
my $x;

#declare array
my @line;
my @id_part;
my @ids;

open (IN,"<","results.clstr") || die "Error cannot read results.fasta file: $!\n";

#####prompt for user input#####
print "Please enter label for this dataset (ex. uclust_097):\n";
$label = <STDIN>;
chomp $label;

$num_otus = qx(grep ">" results.clstr | wc -l);
chomp $num_otus;
print "num_otus: $num_otus\n";#test

$outfile1 = "seed.fasta.list";
open (OUT1,">>",$outfile1) || die "(Error cannot write to .list file: $!\n)";

print OUT1 "$label\t$num_otus";

while (<IN>) {
	$line = $_;
	chomp $line;
	if (/^>/) {
		if (/^>/) {
	#		@line = split (/\|/,$line);
	#		$id_part = $line[2];
#			@id_part = split(/ /,$id_part);
#			$id = $id_part[0];
#			push(@ids,$id);
#			print OUT1 "\t$id";
		}
	}
	else {
		if ($line =~ /\*/) {
			@line = split (/,/,$line);
			$id_part = $line[1];
			$id_part =~ s/\s+>//;
			$id_part =~ s/^>//;#add this
			$id_part =~ s/\.\.\.//;#add this
			@id_part = split(/ /,$id_part);
			$id = $id_part[0];
			#push(@ids,$id);
			#$line =~ />(\w+)\s+/;
			print OUT1 "\t$id";
		}
		else {
			@line = split (/,/,$line);
			$id_part = $line[1];
			$id_part =~ s/\s+>//;
			$id_part =~ s/^>//;#add this
			$id_part =~ s/\.\.\.//;#add this
			@id_part = split(/ /,$id_part);
			$id = $id_part[0];
			#push(@ids,$id);
			print OUT1 ",$id";
		}
	}
}
print OUT1 "\n";
close IN;
close OUT1;

}
