#!/usr/bin/perl
#Nov.12,2010 by Terri Porter
#Script to reformat a fasta file so that it can be properly parsed by agrep.
#Use output to search for 5' primer sequences with/without mismatches allowed
#usage $perl fasta_to_agrep.plx file.fasta primer.txt
#Nov.24,2010 modify to include a remove_newline.plx part
#Nov.25,2010 modify to prompt for input at the start, then automatically run fasta_to_agrep, merge, remove_duplicates, and sort_qual_by_ID_notrim
#Dec.03,2010 modify to automatically run seqtrim in separate directories; rename to primer_sort_clean_pipeline.plx
#my_qsub needs to be run from infoserv head node, change program to accept three arguments at the start instead of prompting for user input
#couldn't get my_qsub to work right so just do below, then wait for finish, then manually set up seqtrim jobs, then process with clean_processing.plx
#USAGE $perl primer_sort.plx *.fna *.qual primers.txt

use strict;
use warnings;

#global var
my $dir;
my $fna_name;
my $qual_name;
my $primer_name;

#global array
my @folder_names;
my @folder_names2;

########################

process_infiles();
print "Starting to remove new line from fna file\n";
remove_new_line_from_fna();
print "Starting to remove empty lines from fna file\n";
remove_empty_lines();
print "Getting raw statistics\n";
get_stats();
print "Starting to convert fasta format for agrep\n";
fasta_to_agrep_format();
print "Starting agrep_search\n";
agrep_search();
print "Starting to merge agrep results for same locus\n";
merge_agrep_results_for_same_locus();
print "Removing possible duplicates from merged files\n";
remove_duplicates_from_merged_files();
print "Removing empty lines from qual file\n";
remove_new_line_from_qual();
print "Removing empty lines from qual file\n";
remove_empty_lines();
print "Sorting phred scores by locus\n";
sort_qual_entries_by_locus();
print "Sorting fastas by qual order\n";
order_fasta_to_qual();
print "Done parsing, now setting up seqtrim jobs\n";
run_seqtrim();

########################

sub process_infiles {

	$dir = qx(pwd);
	chomp $dir;
	$dir =~ s/\/net\/infoserv\/3//;	
	$fna_name = $ARGV[0];
	chomp $fna_name;

	$qual_name = $ARGV[1];
	chomp $qual_name;

	$primer_name = $ARGV[2];
	chomp $primer_name;

}

#########################

sub remove_new_line_from_fna {
	#var
	my $line;

	open (IN,"<",$fna_name) || die ("Error with fna_name: $!\n");
	open (TEMP,">>","temp1.txt") || die ("Error with temp1:$!\n");

	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^>/) {
			print TEMP "\n",$line,"\n";
		}
		else {
			print TEMP $line;
		}
	}

	close IN;
	close TEMP;
}

##########################

sub remove_empty_lines {
	#var
	my $line;

	open (IN,"<","temp1.txt") || die ("Error: $!\n");
	open (TEMP,">>","temp2.txt") || die ("Error:$!\n");

	while (<IN>) {
		$line=$_;
		chomp $line;
		if ($line =~ /\S+/) {
			print TEMP $line,"\n";
		}
	}
	close IN;
	close TEMP;

}

###########################

sub get_stats {
	
	use Statistics::Lite qw(:all);

	#var
	my $line;
	my $flag=0;
	my $seq1;
	my $seq2;
	my $i=0;
	my $length;
	my $min;
	my $max;
	my $mean;
	my $mode;
	my $num;
	my $element;

	#array
	my @seq;
	my @split;
	my @length;

	open (IN, "<","temp2.txt") || die ("Error cannot read from temp1 file: $!\n");
	
	open (OUT,">>","raw.stats") || die ("Error cannot write to raw stats file: $!\n");

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
				$seq1=$seq1.$seq2;
			}
		}
	}
	close IN;
	push (@seq, $seq1);

	while ($seq[$i]){
		@split = split(//,$seq[$i]);
		$length = scalar(@split);
		push (@length,$length);
		@split = ();
		$i++;
	}

	foreach(@length) {
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

##########################
sub fasta_to_agrep_format {
	#var
	my $line;
	my $id;
	my $seq;
	
	open (IN,"<","temp2.txt") || die ("Error opening temp2 file: $!\n");
	open (TEMP,">>","temp3.txt") || die ("Error opening temp3 file: $!\n");
	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^>/) {
			$line =~ /^>(\w{16})\s+/; #####Edit this to work with old 454 IDS#####
			$id = $1;
		}
		elsif ($line =~ /^\w+/) {
			$seq = $line;
			print TEMP "$line|$id\n";
		}
	}
	close IN;
	close TEMP;
}

##########################

sub agrep_search {
	#var
	my $line;
	my $primername;
	my $primersequence;
	my $pattern;
	my $fname;
	my $filename;
	my $i=0;
	my $k=0;
	my $first_letter;
	my $second_letter;
	my $third_letter;
	my $fourth_letter;
	my $fifth_letter;
	my $folder_name;
	my $folder_name_path;

	#array
	my @line;
	my @primername;
	my @primersequence;
	my @output;
	my @fname;
	my @output2;

	open (IN,"<",$primer_name) || die "Error: $!\n";

	while (<IN>) {
		$line = $_;
		chomp $line;
		@line = split(/\t/,$line);
		$primername = $line[0];
		$primersequence = $line[1];
		push(@primername,$primername);
		push(@primersequence,$primersequence);
	}
	close IN;

	while ($primersequence[$i]) {
		$pattern = $primersequence[$i];
		###adjust agrep settings here###
		@output = qx(agrep -i -0 -D1 -I1 "^$pattern" temp3.txt);  #####change -1 to -0 with MIDs#####
		$fname = $primername[$i];
		$filename = $fname.".agrep";
		#make a new folder based on first five letters of primer name
		@fname = split(//,$fname);
		$first_letter = $fname[0];
		$second_letter = $fname[1];
		$third_letter = $fname[2];
		$fourth_letter = $fname[3];
		$fifth_letter = $fname[4];
		$folder_name = $first_letter.$second_letter.$third_letter.$fourth_letter.$fifth_letter;
		mkdir $folder_name;
		push(@folder_names,$folder_name);
		$folder_name_path = $dir."/".$folder_name."/".$filename;

		open (OUT,">>",$folder_name_path) || die "Error: $!\n";
		while ($output[$k]) {
			$line = $output[$k];
			print OUT "$line";
			$k++;
		}
		close OUT;
		$k=0;
		$i++;
	}

	unlink("temp1.txt","temp2.txt","temp3.txt");
	#@output2 = qx(ls -lhrt);
	#print "@output2\n";
}

##########################

sub merge_agrep_results_for_same_locus {

	#var
	my $folder_list;
	my $filename;
	my $file;
	my $pathtofile;
	my $i=0;
	my $locus_name;
	my $locus_name_path;
	my $open_this_dir;

	#array
	my @files;
	my @new_array;
	my @folder_list;
	my @output;

	#print "Enter comma-delimited list of locus folders to analyze (no spaces): ";
	#$folder_list = <STDIN>;
	#chomp $folder_list;
	#@folder_list = split(/,/,$folder_list);
	
	while ($folder_names[$i]) {
		$locus_name = $folder_names[$i];
		$open_this_dir = $dir."/".$locus_name;
		opendir DH, $open_this_dir;
		@files = readdir (DH);
		$locus_name_path = $dir."/".$locus_name.".merged";
		open (OUT, ">>",$locus_name_path) || die ("Cannot open new merged file");

		foreach $file (@files) {
			$filename = $file;
			$pathtofile = $open_this_dir."/".$filename;
			open (FH, $pathtofile) || die ("Cannot open file in array");
			@new_array = <FH>;
			print OUT @new_array;
			close (FH);
		}
		close OUT;
		$i++;
		@new_array=();
		@files=();
	}
	#@output = qx(ls -lhrt);
	#print "@output\n";
}

#########################

sub remove_duplicates_from_merged_files {

	#var
	my $i=0;
	my $line;
	my $seq;
	my $id;
	my $j=0;
	my $k=0;
	my $nonredundant_id;
	my $checkid;
	my $sequence;
	my $flag=0;
	my $filename;
	my $l=0;
	my $filename_new;

	#arrays
	my @file1;
	my @line;
	my @id;
	my @seq;
	my @id_nonredundant;
	my @output;
	my @output2;

	#hash
	my %ID;

	@output = qx(ls | grep '.merged');
	foreach $filename (@output) {
		chomp $filename;
	}
	
	while ($output[$l]) {
		$filename = $output[$l];

		open (IN,"<",$filename) || die "Error opening merged file: $!\n";
		@file1 = <IN>;
		close IN;

		while ($file1[$i]) {
			$line = $file1[$i];
			chomp $line;
			@line = split(/\|/,$line);
			$seq = $line[0];
			$id = $line[1];
			push(@id,$id);
			push(@seq,$seq);
			$i++;
		}

		#dereplicate
		%ID = map {$_,1} @id;
		@id_nonredundant = keys %ID;
		
		$filename_new = $filename.".dereplicated";
		open (OUT, ">>", $filename_new) || die "Error opening dereplicated file: $!\n";
		
		while ($id_nonredundant[$j]) {
			$nonredundant_id = $id_nonredundant[$j];

			while ($id[$k]) {
				$checkid = $id[$k];
				$sequence = $seq[$k];
				if ($nonredundant_id eq $checkid) {
					if ($flag==0) {
						print OUT ">$checkid\n$sequence\n";
						$flag=1;
					}
					elsif ($flag==1) {
						$k++;
						next;
					}
				}
				$k++;
			}
			$flag=0;
			$j++;
			$k=0;
		}
		$l++;
		$i=0;
		$j=0;
		@file1 = ();
		@id = ();
		@seq = ();
		%ID = ();
	}
	@output2 = qx(ls -lhrt);
	print "@output2\n";

}
#########################

sub remove_new_line_from_qual {
	
	#var
	my $line;
	my $id_to_match;
	my $id;
	my $x;
	my $i=1;
	my $trimmed_phred;
	my $flag=0;
	my $w=0;
	my $phred;
	my $header;

	#array
	my @qual;
	my @fna;
	my @ids;
	my @line;

	open (IN,"<",$qual_name) || die ("Error opening qual file: $!\n");
	open (TEMP,">>","temp1.txt") || die ("Error opening temp file: $!\n");

	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^>/) {
			print TEMP "\n",$line,"\n";
			$flag=0;
		}
		elsif ($flag==0) {
			print TEMP $line;
			$flag=1;
		}
		elsif ($flag>0) {
			print TEMP " $line";#need the space here
		}
	}
	close IN;
	close TEMP;
}

###########################

sub sort_qual_entries_by_locus {
	
	#var
	my $filename;
	my $l=0;
	my $line;
	my $id_to_match;
	my $filename_sorted;
	my $w=0;
	my $header;
	my $id;
	my $x;
	my $phred;
	my $basename;

	#array
	my @output;
	my @file;
	my @qual;
	my @ids;
	my @output2;
	my @filename;

	@output = qx(ls | grep '.merged.dereplicated');
	foreach $filename (@output) {
		chomp $filename;
	}

	while ($output[$l]) {
		$filename = $output[$l];

		open (QUAL,"<","temp2.txt") || die ("Error opening temp2 file: $!\n");
		@qual = <QUAL>;

		open (FNA,"<",$filename) || die ("Error opening dereplicated file: $!\n");
		
		while(<FNA>) {
			$line = $_;
			chomp $line;
			if ($line =~ />/) {
				$line =~ />(\w{16})/; #####Edit this to work with old 454 IDS#####
				$id_to_match = $1;
				push(@ids, $id_to_match);
			}
		}
		close FNA;
		
		@filename = split(/\./,$filename);
		$basename = $filename[0];	
		$filename_sorted = $basename.".sorted.qual";
		open (OUT,">>",$filename_sorted) || die ("Error opening sorted file: $!\n");

		while ($qual[$w]) {
			$line = $qual[$w];
			chomp $line;
			if ($line =~ /^>/) {
				$header = $line;
				$header =~ /^>(\w{16})/; #####Edit this to work with old 454 IDS#####
				$id = $1;
				foreach $id_to_match (@ids) {
					if ($id_to_match eq $id) {
						$x = $w+1;
						$phred = $qual[$x];
						chomp $phred;
						print OUT ">$id\n$phred\n";
					}
				}
			}
			$w++
		}
		close QUAL;
		close OUT;
		$l++;
		@qual = ();
		@ids = ();
		$w=0;
	}
	#@output2 = qx(ls -lhrt);
	#print "@output2\n";
	unlink("temp1.txt","temp2.txt");
}

############################

sub order_fasta_to_qual {

	#var
	my $i=0;
	my $line;
	my $id;
	my $j=0;
	my $line2;
	my $id2;
	my $k;
	my $seq;
	my $filename;
	my $l=0;
	my $filename2;
	my $m=0;
	my $basename;
	my $fastaname;
	my $folder_name;
	my $folder_name_path;
	my $qualname;
	my $old_location;
	my $qualname_path;
	my $new_location;
	my $fastaname_sorted;
	my $fastaname_sorted_path;

	#array
	my @fasta;
	my @qual;
	my @output;
	my @output2;
	my @filename;
	my @basenames;

	use File::Copy;

	@output = qx(ls | grep '.merged.dereplicated');
	foreach $filename (@output) {
		chomp $filename;
		@filename = split(/\./,$filename);
		$basename = $filename[0];
		push (@basenames,$basename);
	}

	while ($basenames[$l]) {
		$basename= $basenames[$l];
		$fastaname = $basename.".merged.dereplicated";
		open (FASTA,"<",$fastaname) || die ("Error opening dereplicated fasta: $!\n");
		@fasta = <FASTA>;
		close FASTA;
		
		$folder_name = $basename."_for_seqtrim";
		mkdir $folder_name;
		$folder_name_path = $dir."/".$folder_name;
		push(@folder_names2, $folder_name);		
		$qualname = $basename.".sorted.qual";
		$old_location = $qualname;
		$qualname_path = $folder_name_path."/".$qualname;
		$new_location = $qualname_path;
		move ($old_location, $new_location);

		open (QUAL,"<",$qualname_path) || die ("Error opening sorted qual file: $!\n");
		@qual = <QUAL>;
		close QUAL;
		
		$fastaname_sorted = $basename.".sorted.fasta";
		$fastaname_sorted_path = $folder_name_path."/".$fastaname_sorted;
		open (OUT,">>",$fastaname_sorted_path) || die ("Error opening sorted fasta file: $!\n");

		while ($qual[$i]) {
			$line = $qual[$i];
			chomp $line;
			if ($line =~ /^>/) {
				$line =~ /^>(\w{16})\s*/; #####Edit this to work with old 454 IDS#####
				$id = $1;
				$j=0;
				my $flag=0;
				while ($fasta[$j]) {
					$line2 = $fasta[$j];
					chomp $line2;
					if ($line2 =~ /^>/) {
						$line2 =~ /^>(\w{16})\s*/; #####Edit this to work with old 454 IDS#####
						$id2 = $1;
						if ($id2 eq $id && $flag==0) {
							$k=$j+1;
							$seq = $fasta[$k];
							chomp $seq;
							print OUT ">$id2\n$seq\n";
							$flag=1;
						}
						elsif ($id2 eq $id && $flag==1) {
							$j++;
							next;
						}
					}
					$j++;
				}
			}
			$i++;
		}
		$l++;
		$i=0;
		@fasta = ();
		@qual = ();
	}
	#@output2 = qx(ls -lhrt);
	#print "@output2\n";
}

####################################

sub run_seqtrim {
	
	#var
	my $folder_list;
	my $folder;
	my $fasta_file;
	my $qual_file;
	my $basename;
	my $seqtrim_outfile;
	my $line;
	my $outfile;
	my $fasta_file_path;
	my $qual_file_path;
	my $seqtrim_outfile_path;

	#array
	my @folders;
	my @split;
	my @output;

	#print "Enter comma-separated list of folders that contain sequences for seqtrim read trimming (no spaces)";
	#$folder_list = <STDIN>;
	#@folders = split(/,/,$folder_list);
	
	foreach $folder (@folder_names2) {
		chomp $folder;
		#print $folder."\n";#test
		chdir $folder;
		$fasta_file = qx(ls | grep ".fasta");
		chomp $fasta_file;
		#print $fasta_file."\n";#test
		#$fasta_file_path = $dir."/".$folder."/".$fasta_file;
		$qual_file = qx(ls | grep ".qual");
		chomp $qual_file;
		#print $qual_file."\n";#test
		#$qual_file_path = $dir."/".$folder."/".$qual_file;
		@split = split(/\./,$fasta_file);
		$basename = $split[0];
		$seqtrim_outfile = $basename."_seqtrim.fasta";
		#$seqtrim_outfile_path = $dir."/".$folder."/".$seqtrim_outfile;
		
		system("my_qsub -d 'perl seqtrim.pl -f $fasta_file -q $qual_file -v --sm=80 --QV=20 --QW=10 --arrange vnq -o $seqtrim_outfile'");
		@split=();
		@output=();
		chdir $dir;
	}
}
