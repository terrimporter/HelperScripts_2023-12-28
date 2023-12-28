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
#Jan.18,2010 modified to work with Saina's data
#USAGE $perl primer_sort.plx *.fna *.qual primers.txt mid.txt

#NEWUSAGE $perl Saina_part1_pipeline.plx file.fna file.qual primer.txt mid.txt

use strict;
use warnings;

#global var
my $dir;
my $fna_name;
my $qual_name;
my $primer_name;
my $mid_name;

#global array
my @folder_names_mid;
my @folder_names;
my @folder_names2;
my @folders_to_process;
my @folder_names_primer;

########################

process_infiles();
print "Starting to remove new line from fna and qual files\n";
remove_new_line_from_fna();
print "Starting to remove empty lines from fna and qual files\n";
remove_empty_lines();
print "Getting raw statistics\n";
get_stats();
print "Starting to convert fasta format for agrep\n";
fasta_to_agrep_format();
print "Starting agrep_mid_search and converting back to fasta format\n";
agrep_mid_search();
print "Starting mid trimming and sort .qual file and  remove first 10 phred scores\n";
trim_mid();
print "get exact primer sequences for each library";
get_primer_exact();
print "Convert to fasta to agrep to do search, create specific primer files, Start primer agrep search\n";
agrep_primer_search();
print "Starting to merge agrep results for same locus\n";
merge_agrep_results_for_same_locus();
print "Removing possible duplicates from merged agrep files and converting to fasta file\n";
remove_duplicates_from_merged_files();
#print "Removing empty lines from qual file\n";
#remove_new_line_from_qual();
#print "Removing empty lines from qual file\n";
#remove_empty_lines();
print "Sorting phred scores by locus\n";
sort_qual_entries_by_locus();
print "Sorting fastas by qual order\n";
order_fasta_to_qual();
print "Done parsing, you need to set up manual seqtrim jobs\n";
#run_seqtrim();

########################

sub process_infiles {

	#var
	my $x;
	#array
	my @output;
	
	$dir = qx(pwd);
	chomp $dir;
	$dir =~ s/\/net\/infoserv\/3//;	

	@output = qx(ls);
	foreach $x(@output){
		print $x;
	}
	print "Enter name of raw .fna file:\n";
	
	$fna_name = <STDIN>;
	chomp $fna_name;

	print "Enter name of .qual file:\n";
	$qual_name = <STDIN>;
	chomp $qual_name;

	#$primer_name = $ARGV[2];
	#chomp $primer_name;

	print "Enter name of mid.txt:\n";
	$mid_name = <STDIN>;
	chomp $mid_name;
	
}

#########################

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

##########################

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

##########################

sub agrep_mid_search {
	#var
	my $line;
	my $midname;
	my $midsequence;
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
	my @midname;
	my @midsequence;
	my @output;
	
	open (IN,"<",$mid_name) || die ("Error cannot open mid.txt: $!\n");

	while (<IN>) {
		$line = $_;
		chomp $line;
		@line = split (/\t/,$line);
		$midname = $line[0];
		$midsequence = $line[1];
		push(@midname,$midname);
		push(@midsequence,$midsequence);
	}
	close IN;
	
		
	while ($midsequence[$i]) {
		$pattern = $midsequence[$i];
		###adjust agrep settings here###
		@output = qx(agrep -i -0 -D1 -I1 "^$pattern" temp3.txt);
		$fname=$midname[$i];
		$filename = $fname.".agrep";
		#name a new folder based on MID name
		$folder_name = $fname;
		mkdir $folder_name;
		push(@folder_names_mid,$folder_name);
		$folder_name_path = $dir."/".$folder_name."/".$filename;

		open (OUT,">>",$folder_name_path) || die "Error: $!\n";
	
		while ($output[$k]) {
			$line = $output[$k];
			print OUT "$line";
			$k++;
		}
		close OUT;

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

sub trim_mid {

	#var
	my $x;
	my $folders_to_process;
	my $i=0;
	my $folder_name;
	my $folder_name_path_to_agrep;
	my $line;
	my $outfile;
	my $header;
	my $seq;
	my $j=1;
	my $trimmed_seq;
	my $id_to_match;
	my $id;
	my $y;
	my $match=0;
	my $k=1;
	my $trimmed_phred;
	my $outfile2;

	#array
	my @output;
	my @line;
	my @fna;
	my @ids;

	@output = qx(ls);
	foreach $x (@output) {
		print $x;
	}

	print "Enter names of folders that need to be processed (comma-separated,nospaces):\n";

	$folders_to_process = <STDIN>;
	chomp $folders_to_process;
	@folders_to_process = split(/,/,$folders_to_process);

	while ($folders_to_process[$i]) {
		$folder_name = $folders_to_process[$i];
		$folder_name_path_to_agrep = $dir."/".$folder_name."/".$folder_name.".fasta";
		open (IN, "<", $folder_name_path_to_agrep) || die ("Error: Cannot open file.fasta $!\n");

		while (<IN>) {
			$line = $_;
			chomp $line;
			$outfile = $folder_name_path_to_agrep.".midtrimmed";
			open (OUT,">>", $outfile) || die ("Error: $! Cannot open file.agrep.midtrimmed $!\n");
			if ($line =~ /^>/) {
				$header = $line;
			}
			else {
				$seq = $line;
				@line = split (//,$line);
				while ($j <= 10) {
					shift (@line);
					$j++;
				}
				$trimmed_seq = join('',@line);
				print OUT "$header\n$trimmed_seq\n";
				$j=1;
			}
			$trimmed_seq=();
		}

		close IN;
		close OUT;

		open (QUAL, "<", "temp2.qual") || die ("Error cannot open .qual file: $!\n");
		open (FNA, "<", $outfile) || die ("Error cannot open .midtrimmed file: $!\n");

		while (<FNA>) {
			$line = $_;
			chomp $line;
			if ($line =~ />/) {
				$line =~ />(\w{14})/;
				$id_to_match = $1;
				push(@ids, $id_to_match);
			}
		}
		close FNA;
		
		$outfile2 = $dir."/".$folder_name."/".$folder_name.".qual.midtrimmed";
		open (OUT2,">>",$outfile2) || die ("Error cannot open file.qual: $!\n");

		while (<QUAL>) {
			$line = $_;
			chomp $line;
			if ($line =~ /^>/) {
				$line =~ /^>(\w{14})/;
				$id = $1;
				foreach $y (@ids) {
					if ($y eq $id) {
						$match = 1;
					}
				}
			}
			elsif ($match==1) {
				@line = split(/ /,$line);
				while ($k<=10) {
					shift (@line);
					$k++;
				}
				$trimmed_phred = join(' ',@line);
				print OUT2 ">$id\n$trimmed_phred\n";
				$k=1;
				$match=0;
			}
			$trimmed_phred = ();
		}
		@ids=();
		close QUAL;
		close OUT2;

		$i++;
	}
}

####################

sub get_primer_exact {

	print "\nEnter file name for original primer file (ex. primer.txt):\n";
	$primer_name = <STDIN>;
	chomp $primer_name;

	#declare variable
	my $line;
	my $name;
	my $degenerate_primer;
	my $i=0;
	my $degenerate_base;
	my $new_exact_primer;
	my $old_exact_primer;
	my $x;
	my $exact_base;
	my $y;
	my $z;
	my $scalar;
	my $j=0;

	#declare array
	my @line;
	my @names;
	my @degenerate_primers;
	my @degenerate_bases;
	my @new_exact_primers;
	my @old_exact_primers;

	open (IN,"<",$primer_name) || die ("Error cannot open primer file:$!\n");
	
	while (<IN>) {
		$line = $_;
		chomp $line;
		@line = split(/\t/,$line);
		$name = $line[0];
		$degenerate_primer = $line[1];
		push(@names,$name);
		push(@degenerate_primers,$degenerate_primer);
	}
	close IN;
	
	while ($degenerate_primers[$i]) {
		$degenerate_primer = $degenerate_primers[$i];
		@degenerate_bases = split(//,$degenerate_primer);
		foreach $degenerate_base (@degenerate_bases) {
			if ($degenerate_base =~ /(A|C|T|G)/i) {
				$exact_base = $degenerate_base;
				$scalar = scalar(@new_exact_primers);
				if ($scalar >= 1) {
					foreach $x (@new_exact_primers) {
						push(@old_exact_primers,$x);
					}
					@new_exact_primers = ();
					foreach $old_exact_primer (@old_exact_primers) {
						$new_exact_primer = $old_exact_primer.$exact_base;
						push(@new_exact_primers,$new_exact_primer);
					}
				}
				elsif ($scalar == 0){
					$new_exact_primer = $exact_base;
					push(@new_exact_primers,$new_exact_primer);
				}
				@old_exact_primers = ();
			}
			else {
				foreach $x (@new_exact_primers) {
					push(@old_exact_primers,$x);
				}
				@new_exact_primers=();
				foreach $old_exact_primer (@old_exact_primers) {
					if ($degenerate_base =~ /M/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /R/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /W/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /S/i) {
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base=~ /Y/i) {
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base=~ /K/i) {
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /V/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /H/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /D/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers,$new_exact_primer);
					}
					elsif ($degenerate_base =~ /B/i) {
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers, $new_exact_primer);
					}
					elsif ($degenerate_base =~ /(X|N|I)/i) {
						$new_exact_primer = $old_exact_primer."A";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."C";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."G";
						push(@new_exact_primers,$new_exact_primer);
						$new_exact_primer = $old_exact_primer."T";
						push(@new_exact_primers,$new_exact_primer);
					}
					else {
						print "Unknown nucleotide ambiguity character or other punctuation present in primer sequence.\n";
					}
				}
				@old_exact_primers = ();
			}
		}
		$y = $names[$i];
		my $zz = scalar(@new_exact_primers);
	
		open (PRIMER,">>","primer.exact") || die ("Error cannot open primer.exact: $!\n");
	
		while ($new_exact_primers[$j]) {
			$z = $new_exact_primers[$j];
			$a = $j+1;
			if ($zz>1) {
				my $newprimername = $y."_".$a;
				print PRIMER "$newprimername\t$z\n";				
			}
			else {
				print PRIMER "$y\t$z\n";
			}
			$j++;
		}
		$j=0;
		$i++;
		@new_exact_primers = ();
		close PRIMER;
	}

}

####################
		
sub agrep_primer_search {
	#var
	my $a=0;
	my $folder_name_mid;
	my $folder_name_mid_path;
	my $folder_name_mid_path_agrep;
	my $line;
	my $id;
	my $seq;
	my $primerout;
	my $z;
	my $primername;
	my $primersequence;
	my $mismatches;
	my $i=0;
	my $pattern;	
	my $fname;
	my $filename;
	my $first_letter;
	my $second_letter;
	my $third_letter;
	my $folder_name;
	my $path_for_folder_name;
	my $path_for_filename;
	my $k=0;

	#array
	my @line;
	my @output;
	my @primername;
	my @primersequence;
	my @output2;
	my @fname;
	#my @folder_names_primer;

	while ($folders_to_process[$a]) {
		$folder_name_mid = $folders_to_process[$a];
		print "folder name mid = $folder_name_mid\n";#test
		$folder_name_mid_path = $dir."/".$folder_name_mid."/".$folder_name_mid.".fasta.midtrimmed";
		print "folder name mid path = $folder_name_mid_path\n";#test	
		open (IN,"<",$folder_name_mid_path) || die "Error cannot open file.fasta.midtrimmed:$!\n";				
		#convert to agrep format
		$folder_name_mid_path_agrep = $folder_name_mid_path.".agrep";
		print "folder name mid path agrep = $folder_name_mid_path_agrep\n"; #test
		open (TEMP1,">>",$folder_name_mid_path_agrep) || die ("Error cannot open .temp1 file: $!\n");
			
		while (<IN>) {
			$line = $_;
			chomp $line;
			if ($line =~ /^>/) {
				$line =~ /^>(\w{14})/;
				$id = $1;
			}
			elsif ($line =~ /^\w+/) {
				$seq = $line;
				print TEMP1 "$line|$id\n";
			}
		}
		close IN;
		close TEMP1;

		#create primer file for just this MID
				
		@output = qx(grep "$folder_name_mid" primer.exact);
		$primerout = $dir."/".$folder_name_mid."/"."primer.exact.".$folder_name_mid;
			
		open (OUT,">>",$primerout) || die ("Error cannot open primer.exact.mid:$!\n");
			
		foreach $z (@output) {
			chomp $z;
				print OUT $z."\n";
		}
		close OUT;

		open (IN2,"<",$primerout) || die ("Error cannot open primer.exact.mid2: $!\n");

		while (<IN2>) {
			$line = $_;
			chomp $line;
			@line = split(/\t/,$line);
			$primername = $line[0];
			$primersequence = $line[1];
			push(@primername,$primername);
			push(@primersequence,$primersequence);
		}
		close IN2;

		#print "primer name array: @primername\n";#test
		#print "primer sequence array: @primersequence\n";#test
		
		#do agrep search
		
		print "\nFor library $folder_name_mid, please enter number of allowed agrep mismatches:\n";
		$mismatches = <STDIN>;
		chomp $mismatches;
		
		while ($primersequence[$i]) {
			$pattern = $primersequence[$i];
			###adjust agrep settings here###
			@output2 = qx(agrep -i -$mismatches -D1 -I1 "^$pattern" $folder_name_mid_path_agrep);
			$fname = $primername[$i];
			print "primer fname: $fname\n";#test
			$filename = $fname.".agrep";
			#make a new folder based on primer name
			@fname = split(/_/,$fname);
			$folder_name = $fname[2];
			$path_for_folder_name = $dir."/".$folder_name_mid."/".$folder_name;

			unless (-d $path_for_folder_name) { #only mkdir if it doesn't already exist	
				mkdir $path_for_folder_name;
				push(@folder_names_primer,$path_for_folder_name);
			}
			
			$path_for_filename = $path_for_folder_name."/".$filename;

			open (OUT,">>",$path_for_filename) || die "Error: $!\n";
			while ($output2[$k]) {
				$line = $output2[$k];
				chomp $line;
				print OUT "$line\n";
				$k++;
			}
			close OUT;

			$k=0;
			$i++;
		}
		$i=0;
		@primername=();
		@primersequence=();
		$a++;
	}
	print "#####\narray folder names primer: @folder_names_primer\n";#test
}

##########################

sub merge_agrep_results_for_same_locus {

	#var
	my $current_folder;
	my $y;
	my $folder;
	my $outfile_name;
	my $l=0;
	my $j=0;
	my $agrep_file;

	#array
	my @output3;
	my @agrep_files;
	my @current_folder;
	my @new_array;
	
	while ($folder_names_primer[$l]) {
		$current_folder = $folder_names_primer[$l];
		chdir ("$current_folder");	
		@output3 = qx(ls | grep ".agrep");
		foreach $y (@output3) {
			chomp $y;
			push (@agrep_files,$y);
		}
		
		print "array agrep files: @agrep_files\n";#test
		
		while ($agrep_files[$j]) {
			@current_folder = split(/\//,$current_folder);
			$folder = pop(@current_folder);
			$outfile_name = $folder.".merged";
			
			open (OUT,">>",$outfile_name) || die ("Error cannot open .merged file: $!\n");
		
			$agrep_file = $agrep_files[$j];
		
			open (FH, "<", $agrep_file) || die ("Cannot open agrep file: $!\n");
			
			@new_array = <FH>;
			print OUT @new_array;
			close (FH);
			$j++;
		}
		$j=0;
		close OUT;
		@new_array=();
		$l++;
		@agrep_files=();
	}
	$l=0;
	chdir ("/home/terri/Saina/");
}

#########################

sub remove_duplicates_from_merged_files {

	#var
	my $m=0;
	my $current_folder;
	my $output4;
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
	
	print "folder names primer: @folder_names_primer\n";#test

	while ($folder_names_primer[$m]) {
		$current_folder = $folder_names_primer[$m];
		chdir ("$current_folder");
		#my $test = qx(pwd);#test
		#print "###\npwd: $test\n###\n";#test
		$output4 = qx(ls | grep ".merged");
		chomp $output4;

		open (IN,"<",$output4) || die "Error opening merged file: $!\n";
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
		
		my @output4 = split(/\./,$output4);
		$filename_new = $output4[0].".dereplicated";
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
		$m++;
		$i=0;
		$j=0;
		@file1 = ();
		@id = ();
		@seq = ();
		%ID = ();
	}
	chdir ("/home/terri/Saina/");
	#@output2 = qx(ls -lhrt);
	#print "@output2\n";

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

	while ($folder_names_primer[$n]) {
		$current_folder = $folder_names_primer[$n];
		chdir ("$current_folder");
		$output5 = qx(ls | grep ".dereplicated");
		chomp $output5;
		$filename = $output5;

		open (FNA,"<",$filename) || die ("Error opening dereplicated file: $!\n");
		@fna = <FNA>;
		close FNA;
		#print "fna array: @fna\n";#test	
		@current_folder = split(/\//,$current_folder);
		my $scalar = scalar(@current_folder);
		my $splice = $scalar-1;
		my @up_dir = splice(@current_folder,0,$splice);
		$up_dir = join('/',@up_dir);
		my $up_dir2 = $up_dir."/";
		chdir ("$up_dir2");
		#my $test2 = qx(pwd);
		#print "pwd updir : $test2\n";
		$output6 = qx(ls | grep ".qual.midtrimmed");
		chomp $output6; 
		open (QUAL,"<",$output6) || die ("Error opening .qual.midtrimmed file: $!\n");
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
	}
	chdir("/home/terri/Saina/");
}

############################

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

	use File::Copy;

	while ($folder_names_primer[$p]) {
		$current_folder = $folder_names_primer[$p];
		chdir ("$current_folder");
		$output7 = qx(ls | grep ".dereplicated");
		chomp $output7;
		
		open (FASTA,"<",$output7) || die ("Error opening dereplicated fasta: $!\n");
		@fasta = <FASTA>;
		close FASTA;
		
		@output7 = split(/\./,$output7);
		$basename = $output7[0];
		$folder_name = $basename."_for_seqtrim";
		mkdir $folder_name;
		$output8 = qx(ls | grep ".sorted.qual");
		chomp $output8;
		
		open (QUAL,"<",$output8) || die ("Error opening sorted qual file: $!\n");
		@qual = <QUAL>;
		close QUAL;
		
		$folder_name_path = $current_folder."/".$folder_name;
		push(@folder_names2, $folder_name);		
		$old_location = $current_folder."/".$output8;
		$new_location = $folder_name_path."/".$output8;
		copy ($old_location, $new_location);
		
		$fastaname_sorted = $basename.".sorted.fasta";
		$fastaname_sorted_path = $folder_name_path."/".$fastaname_sorted;
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
	}
	chdir("/home/terri/Saina/");
	unlink("temp1.txt","temp.qual","temp2.txt","temp2.qual","temp3.txt");
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
		system("exit");		
		system("my_qsub -d 'perl seqtrim.pl -f $fasta_file -q $qual_file -v --sm=80 --QV=20 --QW=10 --arrange vnq -o $seqtrim_outfile'");
		@split=();
		@output=();
		chdir $dir;
	}
}
