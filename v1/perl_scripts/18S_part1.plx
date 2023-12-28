#!/usr/bin/perl
#Script to process Lepidoptera 18S files, clean up raw, get stats for raw
#Get list of primers, create a workfolder for each primer
#convert raw sequence file to one that can be parsed by agrep
#parse .agrep by each primer and put results into appropriate folder, and turn back into fasta formatted file
#sort phred scores by locus
#order fasta to qual
#manually set up seqtrim jobs
#usage $perl 18S_part1.plx

use strict;
use warnings;

#global var
my $dir;
my $fna_name;
my $qual_name;
my $primer_name;

#global array
my @folder_names;

##########

process_infiles();

print "Starting to remove new line from fna and qual files\n";
remove_new_line_from_fna();

print "Starting to remove empty lines from fna and qual files\n";
remove_empty_lines();

print "Getting raw statistics\n";
get_stats();

print "Starting to convert fasta to format searchable by agrep\n";
fasta_to_agrep_format();

print "Create specific primer files, start primer agrep search\n";
agrep_primer_search();

print "Sorting phred scores by primer\n";
sort_qual_entries_by_primer();

print "Sort fastas by qual order\n";
order_fasta_to_qual();

print "Done parsing, set up seqtrim jobs manually.\n";

####################

sub process_infiles {

	#var
	my $line;
	
	#array
	my @output;
	
	$dir = qx(pwd);
	chomp $dir;
	$dir =~ s/\/net\/infoserv\/3//;	

	#display contents of directory
	@output = qx(ls);
	foreach $line (@output){
		print $line;
	}

	print "Enter name of raw .fna file:\n";
	
	$fna_name = <STDIN>;
	chomp $fna_name;

	print "Enter name of .qual file:\n";
	$qual_name = <STDIN>;
	chomp $qual_name;

	print "Enter name of primer file:\n";
	$primer_name = <STDIN>;
	chomp $primer_name;

}

####################

sub remove_new_line_from_fna {

	#var
	my $line;
	my $flag=0;
	
	#remove newline from raw.fna
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
	
	#remove newline from raw.qual
	open (IN2,"<",$qual_name) || die ("Error with qual_name: $!\n");
	open (TEMP2,">>","temp.qual") || die ("Error with temp.qual: $!\n");

	while (<IN2>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^>/) {
			print TEMP2 "\n",$line,"\n";
			$flag=0;
		}
		elsif ($flag==0) {
			print TEMP2 $line;
			$flag=1;
		}
		elsif ($flag>0) {
			print TEMP2 " $line";
		}
	}
	close IN2;
	close TEMP2;
}

####################

sub remove_empty_lines {
	#var
	my $line;
	
	#remove first empty line from temp1.txt
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

	#remove first empty line from temp.qual
	
	open (IN2,"<","temp.qual") || die ("Error: $!\n");
	open (TEMP2,">>","temp2.qual") || die ("Error: $!\n");

	while (<IN2>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^\S+/) {
			print TEMP2 $line,"\n";
		}
	}
	close IN2;
	close TEMP2;

}

####################

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

#####################

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
			$line =~ /^>(\w{14})\s+/;
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

####################

sub agrep_primer_search {

	#var
	my $line;
	my $primername;
	my $primersequence;
	my $i=0;
	my $pattern;
	my $fname;
	my $filename;
	my $folder_name;
	my $folder_name_path;
	my $k=0;
	my $folder_name_path2;
	my $seq;
	my $id;

	#array
	my @line;
	my @primername;
	my @primersequence;
	my @output;
	
	open (IN,"<",$primer_name) || die ("Error cannot open primer file: $!\n");

	while (<IN>) {
		$line = $_;
		chomp $line;
		@line = split (/\t/,$line);
		$primername = $line[0];
		$primersequence = $line[1];
		push(@primername,$primername);
		push(@primersequence,$primersequence);
	}
	close IN;
	
		
	while ($primersequence[$i]) {
		$pattern = $primersequence[$i];

		###adjust agrep settings here###
		@output = qx(agrep -i -0 -D1 -I1 "^$pattern" temp3.txt);
		$fname=$primername[$i];
		$filename = $fname.".agrep";
		#name a new folder based on primer name
		$folder_name = $fname;
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
		
		#reformat agrep output into a fasta formated file
		open (IN2,"<",$folder_name_path) || die ("Error: $!\n");
		$folder_name_path2 = $dir."/".$folder_name."/".$fname.".fasta";
		open (OUT2, ">>", $folder_name_path2) || die ("Error: $!\n");

		while (<IN2>) {
			$line = $_;
			chomp $line;
			@line = split (/\|/,$line);
			$seq = $line[0];
			$id = $line[1];
			print OUT2 ">$id\n$seq\n";
		}
		close IN2;
		close OUT2;
		$k=0;
		$i++;
	}
}

####################

sub sort_qual_entries_by_primer {
	
	#var
	my $n=0;
	my $current_folder;
	my $output5;
	my $filename;
	my $remove_dir;
	my $remove_dir2;
	my $up_dir;
	my $output6;
	my $o=0;
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
	my @fna;
	my @current_folder;
	my @output;
	my @file;
	my @qual;
	my @ids;
	my @output2;
	my @filename;

	while ($folder_names[$n]) {
		$current_folder = $folder_names[$n];
		print "###\ncurrent folder: $current_folder\n";#test
		chdir ("$current_folder");
		$output5 = qx(ls | grep ".fasta");
		chomp $output5;
		$filename = $output5;
		print "###\noutput5 filename: $filename\n";#test

		open (FNA,"<",$filename) || die ("Error opening fasta file: $!\n");
		@fna = <FNA>;
		close FNA;
	
		chdir ("/home/terri/Lepidoptera/18S/"); 

		open (QUAL,"<","temp2.qual") || die ("Error opening temp2.qual file: $!\n");
		@qual = <QUAL>;
		close QUAL;
		
		while($fna[$o]) {
			$line = $fna[$o];
			chomp $line;
			if ($line =~ />/) {
				$line =~ />(\w{14})/;
				$id_to_match = $1;
				push(@ids, $id_to_match);
			}
			$o++;
		}
		$o=0;
		
		@filename = split(/\./,$filename);
		$basename = $filename[0];	
		$filename_sorted = $basename.".sorted.qual";
		chdir ("$current_folder");
		open (OUT,">>",$filename_sorted) || die ("Error opening sorted file: $!\n");

		while ($qual[$w]) {
			$line = $qual[$w];
			chomp $line;
			if ($line =~ /^>/) {
				$header = $line;
				$header =~ /^>(\w{14})/;
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
		close OUT;
		$l++;
		@qual = ();
		@fna = ();
		@ids = ();
		$w=0;
		$n++;
		chdir("/home/terri/Lepidoptera/18S/");
	}
	chdir("/home/terri/Lepidoptera/18S/");
}

####################

sub order_fasta_to_qual {

	#var
	my $p=0;
	my $current_folder;
	my $output7;
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
	my $output8;
	my $folder_name_path;
	my $qualname;
	my $old_location;
	my $qualname_path;
	my $new_location;
	my $fastaname_sorted;
	my $fastaname_sorted_path;

	#array
	my @fasta;
	my @output7;
	my @qual;
	my @output;
	my @output2;
	my @filename;
	my @basenames;
	my @folder_names2;

	use File::Copy;

	while ($folder_names[$p]) {
		$current_folder = $folder_names[$p];
		print "###\ncurrent folder: $current_folder\n";#test
		chdir ("$current_folder");
		$output7 = qx(ls | grep ".fasta");
		chomp $output7;
		print "###\noutput7: $output7\n";#test
		
		open (FASTA,"<",$output7) || die ("Error opening file.fasta: $!\n");
		@fasta = <FASTA>;
		close FASTA;
		
		@output7 = split(/\./,$output7);
		$basename = $output7[0];
		print "###\nbasename: $basename\n";#test
		$folder_name = $basename."_for_seqtrim";
		mkdir $folder_name;
		$output8 = qx(ls | grep ".sorted.qual");
		chomp $output8;
		
		open (QUAL,"<",$output8) || die ("Error opening sorted qual file: $!\n");
		@qual = <QUAL>;
		close QUAL;
		
		$folder_name_path = $folder_name;
		push(@folder_names2, $folder_name);		
		$old_location = $output8;
		print "###\nold location: $old_location\n";#test
		$new_location = $folder_name_path."/".$output8;
		print "new location: $new_location\n";#test
		copy ($old_location, $new_location);
		
		$fastaname_sorted = $basename.".sorted.fasta";
		$fastaname_sorted_path = $folder_name."/".$fastaname_sorted;
		print "###\nfastaname sorted path: $fastaname_sorted_path\n";#test
		open (OUT,">>",$fastaname_sorted_path) || die ("Error opening sorted fasta file: $!\n");

		while ($qual[$i]) {
			$line = $qual[$i];
			chomp $line;
			if ($line =~ /^>/) {
				$line =~ /^>(\w{14})\s*/;
				$id = $1;
				$j=0;
				my $flag=0;
				while ($fasta[$j]) {
					$line2 = $fasta[$j];
					chomp $line2;
					if ($line2 =~ /^>/) {
						$line2 =~ /^>(\w{14})\s*/;
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
		$p++;
		chdir("/home/terri/Lepidoptera/18S/");
	}
	chdir("/home/terri/Lepidoptera/18S/");
	unlink("temp1.txt","temp.qual","temp2.txt","temp2.qual","temp3.txt");
}

#####################

