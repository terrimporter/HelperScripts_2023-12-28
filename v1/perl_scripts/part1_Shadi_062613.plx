#!/usr/bin/perl
##May 27, 2013 edited to process Shadi's files specifically
##### Assumes primers.txt is in dir/.  First five letters of primer name indicates folder name (and potential pooling).  Further letters can distinguish among primers but not used for sorting into folders.
##### Current working directory is hardcoded for now.  Edit as needed.
#
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
use File::Copy;

#global vars
my $primer_file = 'primers.txt'; ### hard coded, update this as needed
my $pwd = '/1/scratch/terri/Biomonitoring_soil_062513/';  ### hard coded, update this as needed
my $grep_term;
my $file_ext;
my $file_extb;
my $search_term;
my $file_ext2;
my $output;
my $search_term2;
my $file_ext3;

#global arrays
my @folder_names;
my @output;

#global hash
my %hash; #indexed by id, for all ids from every qual file

########################

print "Fix filenames, no spaces, no dashes, sites only.\n";
fix_filenames();

print "Use MOTHUR v1.25.1 to extract trimmed sequences and qualities from sff files.\n";
extract_sff_files();

print "Removing new lines from sequence part of fasta sequence files.\n";
$grep_term = "fasta";
$file_ext = ".fasta.newlineremoved";
grep_dir1(\$grep_term,\$file_ext);

print "Starting to remove empty lines from fasta file.\n";
$grep_term = "newlineremoved";
$file_ext = ".fasta.newlineremoved.emptyremoved";
grep_dir2(\$grep_term,\$file_ext);

$output = qx(mkdir fasta_newlineremoved);
$output = qx(find *.newlineremoved -maxdepth 1 -type f -exec mv {} fasta_newlineremoved \\;);

print "Getting raw statistics.\n";
$grep_term = "emptyremoved";
$file_ext = ".fasta.newlineremoved.emptyremoved.stats";
grep_dir3(\$grep_term,\$file_ext);

$output = qx(mkdir raw_stats);
$output = qx(find *.stats -maxdepth 1 -type f -exec mv {} raw_stats \\;);

print "Starting to convert fasta format for agrep.\n";
$grep_term = "emptyremoved";
$file_ext = ".agrep";
grep_dir4(\$grep_term,\$file_ext);

print "Starting agrep_search.\n";
$grep_term = "agrep";
$file_ext = ".agrep.search";
grep_dir5(\$grep_term,\$file_ext);

$output = qx(mkdir for_agrep);
$output = qx(find *.agrep -maxdepth 1 -type f -exec mv {} for_agrep \\;);
$output = qx(mkdir fasta_emptyremoved);
$output = qx(find *.emptyremoved -maxdepth 1 -type f -exec mv {} fasta_emptyremoved \\;);
$output = qx(mkdir extracted_fasta);
$output = qx(find *.fasta -maxdepth 1 -type f -exec mv {} extracted_fasta \\;);

print "Convert agrep back to fasta format.\n";
$search_term = "search\$";
$file_ext = ".fasta";
process_by_marker(\$search_term,\$file_ext);

print "Remove newline from sequence part of qual files.\n";
$grep_term = "qual";
$file_ext = ".qual.newlineremoved";
grep_dir1(\$grep_term,\$file_ext);

print "Starting to remove empty lines from qual files.\n";
$grep_term = "newlineremoved";
$file_ext = ".qual.newlineremoved.emptyremoved";
grep_dir2(\$grep_term,\$file_ext);

$output = qx(mkdir qual_newlineremoved);
$output = qx(find *.newlineremoved -maxdepth 1 -type f -exec mv {} qual_newlineremoved \\;);
$output = qx(mkdir extracted_qual);
$output = qx(find *.qual -maxdepth 1 -type f -exec mv {} extracted_qual \\;);

print "Grab associated qual for each fasta.\n";
$search_term = "fasta\$";
$file_ext = ".qual.newlineremoved.emptyremoved";
$file_ext2 = ".qual";
process_by_marker2(\$search_term,\$file_ext,\$file_ext2);

$output = qx(mkdir qual_emptyremoved);
$output = qx(find *.qual.newlineremoved.emptyremoved -maxdepth 1 -type f -exec mv {} qual_emptyremoved \\;);

print "Use MOTHUR v1.25.1 to create fastq files.\n";
$search_term = "fasta\$";
process_by_marker3(\$search_term);

print "Use CUTADAPT v1.1 to remove primers and quality trim; Re-create fasta files.\n";
$search_term = "fastq\$";
$file_ext = ".1.fastq";
$file_extb = ".fastq.trimmed";
$file_ext2 = ".fasta.trimmed";
$file_ext3 = ".fasta.trimmed.stats";
process_by_marker4(\$search_term,\$file_ext,\$file_extb,\$file_ext2,\$file_ext3);

$output = qx(mkdir mothur_logfiles);
$output = qx(find mothur.* -maxdepth 1 -type f -exec mv {} mothur_logfiles \\;);

print "Move files into directories for each marker.\n";
process_by_marker5();

print "Part 1 of the pipeline complete, use part 2 to proceed.\n";

####################

sub grep_dir1 {

#var
$grep_term = ${$_[0]};
$file_ext = ${$_[1]};
my $i=0;
my $file;
my $basename;
my $outfile;

#array
my @output;
my @file;

@output = qx(ls | grep $grep_term);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$basename = $file[0];
	$outfile = $basename.$file_ext;
	remove_new_line_from_fasta(\$file,\$outfile);
	$i++;
}
$i=0;

}

####################

sub grep_dir2 {

#var
$grep_term = ${$_[0]};
$file_ext = ${$_[1]};
my $i=0;
my $file;
my $basename;
my $outfile;

#array
my @output;
my @file;

@output = qx(ls | grep $grep_term);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$basename = $file[0];
	$outfile = $basename.$file_ext;
	remove_empty_lines(\$file,\$outfile);
	$i++;
}
$i=0;

}

####################

sub grep_dir3 {

#var
$grep_term = ${$_[0]};
$file_ext = ${$_[1]};
my $i=0;
my $file;
my $basename;
my $outfile;

#array
my @output;
my @file;

@output = qx(ls | grep $grep_term);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$basename = $file[0];
	$outfile = $basename.$file_ext;
	get_stats(\$file,\$outfile);
	$i++;
}
$i=0;

}

####################

sub grep_dir4 {

#var
$grep_term = ${$_[0]};
$file_ext = ${$_[1]};
my $i=0;
my $file;
my $basename;
my $outfile;

#array
my @output;
my @file;

@output = qx(ls | grep $grep_term);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$basename = $file[0];
	$outfile = $basename.$file_ext;
	fasta_to_agrep_format(\$file,\$outfile);
	$i++;
}
$i=0;

}

####################

sub grep_dir5 {

#var
$grep_term = ${$_[0]};
$file_ext = ${$_[1]};
my $i=0;
my $file;
my $basename;
my $outfile;

#array
my @output;
my @file;

@output = qx(ls | grep $grep_term);

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;
	@file = split(/\./,$file);
	$basename = $file[0];
	$outfile = $basename.$file_ext;
	agrep_search(\$file,\$outfile);
	$i++;
}
$i=0;

}

####################

sub process_by_marker {

#var
$search_term = ${$_[0]};
$file_ext = ${$_[1]};
my $i=0;
my $dir;
my $path_dir;
my $scalar;
my $k=0;
my $file;
my $basename;
my $path_file;
my $outfile; #fasta file
my $j=0;
my $line;
my $seq;
my $id;

#array
my @files; #list of files in dir
my @file; #filename split
my @in; #agrep file
my @line;

while ($folder_names[$i]) {
		$dir = $folder_names[$i];
		chomp $dir;
		print "Processing folder: $dir ... \n";

	if ($dir !~ /^\./) {
		$path_dir = $pwd.$dir."/";
		opendir(DIR,$path_dir) || die "Cannot open dir $path_dir:$!\n";
		@files = readdir(DIR);
		$scalar = scalar(@files);
#		print "MarkerFolders:@files\n";

		while ($files[$k]) {
			$file = $files[$k];
			chomp $file;

			if ($file !~ /^\./ && $file =~ /$search_term/) {
				@file = split(/\./,$file);
				$basename = $file[0];

				$path_file = $path_dir.$file;
				open(FILE, "<", $path_file) || die "Error cannot open agrep file in primer directory: $!\n";
				@in = <FILE>;
				close FILE;
		
				$outfile = $path_dir.$basename.$file_ext;
				open (OUT, ">>", $outfile) || die "Error cannot open fasta outfile in primer directory: $!\n";

				while ($in[$j]) {
					$line = $in[$j];
					chomp $line;

					@line = split(/\|/,$line);
					$seq = $line[0];
					$id = $line[1];
					print OUT ">$id\n$seq\n";

					$j++;
					$line=();
					@line=();
					$id=();
					$seq=();
				}
				$j=0;
				close OUT;
			}
			$k++;
			@file=();
			$basename=();
			@in=();
		}
		$k=0;
	}
	@files=();
	$dir=();
	$i++;
}
$i=0;

}

####################

sub process_by_marker2 {

#var
$search_term = ${$_[0]};
$file_ext = ${$_[1]};
$file_ext2 = ${$_[2]};

my $i=0;
my $dir;
my $path_dir;
my $k=0;
my $file;
my $basename;
my $path_file;
my $path_file2;
my $l=0;
my $j=0;
my $line;
my $header;
my $id;
my $m;
my $nextline;
my $phred;
my $outfile; #qual file

#array
my @files; #list of files in dir
my @file; #filename split
my @in; #fasta file
my @in2; #qual file
my @line;

#hash
my %hash; #indexed by id, for all ids from every qual file

while ($folder_names[$i]) {
		$dir = $folder_names[$i];
		chomp $dir;

	if ($dir !~ /^\./) {
		$path_dir = $pwd.$dir."/";
		opendir(DIR,$path_dir) || die "Cannot open dir $path_dir:$!\n";
		@files = readdir(DIR);

		while ($files[$k]) {
			$file = $files[$k];
			chomp $file;

			if ($file !~ /^\./ && $file =~ /$search_term/) {
				@file = split(/\./,$file);
				$basename = $file[0];

				$path_file = $path_dir.$file;
				open(FILE, "<", $path_file) || die "Error cannot open fasta file in primer directory: $!\n";
				@in = <FILE>;
				close FILE;
		
				$path_file2 = $pwd.$basename.$file_ext;
				open(QUAL, "<", $path_file2) || die "Error cannot open qual file: $!\n";
				@in2 = <QUAL>;
				close QUAL;
		
				#hash qual file
				while ($in2[$l]) {
					$line = $in2[$l];
					chomp $line;

					if ($line =~ /^>/) {
						$header = $line;
						$header =~ />(\w{14,16})/;
						$id = $1;
						$m = $l+1;
						$nextline = $in2[$m];
						chomp $nextline;
						$phred = $nextline;
						$hash{$id} = $phred;
					}
					$l++;
					$line=();
					$id=();
					$m=();
					$nextline=();
					$phred=();
				}
				$l=0;

				$outfile = $path_dir.$basename.$file_ext2;
				open (OUT, ">>", $outfile) || die "Error can't open qual outfile in primer directory: $!\n";

				while ($in[$j]) {
					$line = $in[$j];
					chomp $line;

					if ($line =~ /^>/) {
						$header = $line;
						$header =~ /^>(\w{14,16})/;
						$id = $1;
				
						if (exists $hash{$id}) {
							$phred = $hash{$id};
							print OUT ">$id\n$phred\n";
						}
					}
					$j++;
					$line=();
					$header=();
					$id=();
					$phred=();
				}
				$j=0;
				close OUT;
			}
			$k++;
			@file=();
			$basename=();
			@in=();
		}
		$k=0;
		$i++;
		$dir=();
	}
}
$i=0;

}

####################

sub process_by_marker3 {

#var
$search_term = ${$_[0]};
my $fastafile;
my $qualfile;

my $i=0;
my $dir;
my $path_dir;
my $k=0;
my $file;
my $basename;
my $output;

#array
my @files; #list of files in dir
my @file; #filename split

while ($folder_names[$i]) {
	$dir = $folder_names[$i];
	chomp $dir;

	if ($dir !~ /^\./) {
		$path_dir = $pwd.$dir."/";
		opendir(DIR,$path_dir) || die "Cannot open dir $path_dir:$!\n";
		@files = readdir(DIR);

		while ($files[$k]) {
			$file = $files[$k];
			chomp $file;
#			print "file:$file\n";

			if ($file !~ /^\./ && $file =~ /$search_term/) {
#				print "file:$file\n";
				@file = split(/\./,$file);
				$basename = $file[0];
				$fastafile = $path_dir.$file;	
				$qualfile = $path_dir.$basename.".qual";	
				$output = qx(mothur "#make.fastq(fasta=$fastafile,qfile=$qualfile); quit()");
			}
			$k++;
			$file=();
			@file=();
			$basename=();
			$fastafile=();
			$qualfile=();
			$output=();
		}
		$k=0;
		$path_dir=();
		@files=();
	}
	$i++;
	$dir=();
}
$i=0;

}

####################

sub process_by_marker4 {

#var
$search_term = ${$_[0]};
$file_ext = ${$_[1]}; #first cutadapt output file
$file_extb = ${$_[2]}; #second cutadapt output file
$file_ext2 = ${$_[3]};
$file_ext3 = ${$_[4]};
my $dir;
my $output;
my $primer;
my $infile;
my $outfile;
my $i=0;
my $path_dir;
my $k=0;
my $file;
my $basename;
my $fastafile;
my $j=0;
my $line;
my $l;
my $nextline;
my $outfile3;
my $string;
my $offset;
my $fragment;
my $primer2;
my $infile2;
my $outfile2;

#array
my @output;
my @files; #list of files in dir
my @file; #filename split
my @in;

while ($folder_names[$i]) {
	$dir = $folder_names[$i];
	chomp $dir; 

	#get primer that matches dir name
	$output = qx(grep $dir $primer_file);
	@output = split(/\t/,$output);
	$primer = $output[1];
	chomp $primer;

	#get the second primer for this marker too
	$string = $dir;
	$offset = length($string)-1;

	if ($string =~ /F$/) {
		$fragment = substr $string, $offset, 1, 'R';
		print "dir:$dir\tnew string:$string\n";
	}
	elsif ($string =~ /R$/) {
		$fragment = substr $string, $offset, 1, 'F';
		print "dir:$dir\tnew string:$string\n";
	}

	$output = qx(grep $string $primer_file);
	@output = split(/\t/,$output);
	$primer2 = $output[1];
	chomp $primer2;

	if ($dir !~ /^\./) {
		$path_dir = $pwd.$dir."/";
		opendir(DIR,$path_dir) || die "Cannot open dir $path_dir:$!\n";
		@files = readdir(DIR);

		while ($files[$k]) {
			$file = $files[$k];
			chomp $file;
#			print "file:$file\n";
			
			if ($file !~ /^\./ && $file =~ /$search_term/) {
				@file = split(/\./,$file);
				$basename = $file[0];
				$infile = $path_dir.$file;
				$outfile = $path_dir.$basename.$file_ext;
				
				$output = qx(cutadapt -b $primer -e 0.1 -O 15 -n 2 $infile > $outfile);
				# -b infile, check 3' and 5' end for a match, best one gets removed
				# -o outfile
				# -e max error rate (num error / read length) allowed for adapter match
				# ### I allowed up to 2 errors during the sorting by primer up to this point
				# ### Need to allow for 'K' ambiguity in 16V6F primer (20bp, 2 errors allowed with -e set to 0.1)
				# ### 0.1 (default) will allow for 1.7 to 2.6 errors depending on the length of the primer
				# -m min read length
				# -M max read length
				# -q quality cutoff 
				# ### same as BWA, substract CUTOFF from qualities; 
				# ### compute partial sums from all indices to the end of seq; 
				# ### cut seq at index where sum is minimal
				# -O minimum overlap length required for adapter (primer) removal (shortest primer 17 bp)
				# -n try to remove adapters at most this many times
				
				#need to check for fwd and rev primer for each marker!
				$infile2 = $outfile;
				$outfile2 = $path_dir.$basename.$file_extb;

				$output = qx(cutadapt -b $primer2 -e 0.1 -m 100 -M 400 -q 20 -O 15 -n 2 $infile2 > $outfile2);
				
				#extract fasta and qual from fastq
				print "Converting fastq back to fasta and qual trimmed files\n";
				open (IN, "<", $outfile2) || die "Error cannot open outfile from cutadapt: $!\n";
				@in = <IN>;
				close IN;
				
				$fastafile = $path_dir.$basename.$file_ext2;
				open (FASTA, ">>", $fastafile ) || die "Error cannot open trimmed fasta file: $!\n";
		
				while($in[$j]) {
					$line = $in[$j];
					chomp $line;

					if ($line =~ /^\@{1}/) {
						$line =~ s/^\@/\>/;
						print FASTA $line."\n";
						$l = $j+1;
						$nextline = $in[$l];
						chomp $nextline;
						print FASTA $nextline."\n";
						$j+=4;
					}
					$line=();
					$l=();
					$nextline=();
				}
				$j=0;
				close FASTA;
				@in=();
				$outfile3 = $path_dir.$basename.$file_ext3;
				print "Getting trimmed read stats.\n";
				get_stats(\$fastafile,\$outfile3);
				$fastafile=();
				
			}
			$k++;
			$file=();
			@file=();
			$basename=();
			$infile=();
			$outfile=();
			$output=();
		}
		$k=0;
		$path_dir=();
		@files=();
	}
	$i++;
	$dir=();
}
$i=0;

}

####################

sub process_by_marker5 {

#var
my $dir;
my $output;
my $i=0;
my $path_dir;

#array
my @output;

while ($folder_names[$i]) {
	$dir = $folder_names[$i];
	chomp $dir;

	if ($dir !~ /^\./) {
		$path_dir = $pwd.$dir."/";
			
		my $testdir = $path_dir."stats"; #check for one only
		unless (-d $testdir) {
			$output = qx(mkdir {$path_dir/stats,$path_dir/trimmed,$path_dir/fastq,$path_dir/qual,$path_dir/fasta,$path_dir/from_agrep});
		}
		my $search = $path_dir."*.stats";
		my $target = $path_dir."stats";
		$output = qx(find $search -maxdepth 1 -type f -exec mv {} $target \\;);

		$search = $path_dir."*.trimmed";
		$target = $path_dir."trimmed";
		$output = qx(find $search -maxdepth 1 -type f -exec mv {} $target \\;);
			
		$search = $path_dir."*.fastq";
		$target = $path_dir."fastq";
		$output = qx(find $search -maxdepth 1 -type f -exec mv {} $target \\;);
			
		$search = $path_dir."*.qual";
		$target = $path_dir."qual";
		$output = qx(find $search -maxdepth 1 -type f -exec mv {} $target \\;);
			
		$search = $path_dir."*.fasta";
		$target = $path_dir."fasta";
		$output = qx(find $search -maxdepth 1 -type f -exec mv {} $target \\;);
			
		$search = $path_dir."*.search";
		$target = $path_dir."from_agrep";
		$output = qx(find $search -maxdepth 1 -type f -exec mv {} $target \\;);
				
		$path_dir=();
	}
	$i++;
	$dir=();
	$output=();
}
$i=0;

}

####################

sub fix_filenames {

#declare var
my $i=0;
my $name_original;
my $name_quotes;
my $name_new;
my $output;

#declare arrays
my @filenames;

@filenames = qx(ls | grep sff);

while ($filenames[$i]) {
	$name_original = $filenames[$i];
	chomp $name_original;

	if ($name_original !~ /^\.+/) {
#		$name_quotes = "'".$name_original."'";
		$name_new = $name_original;

		if ($name_new =~ /(PAD|WC)/) { 
			$name_new =~ s/\s+//g; #remove spaces
			$name_new =~ s/\-//g; #remove dashes
			$name_new =~ /\S*((PAD|WC)\d+\S+)/;
			$name_new = $1;
		}
		
		move($name_original, $name_new) || die (qq{failed to move $name_original -> $name_new});
	}

	$i++;
	$name_original=();
	$name_quotes=();
	$name_new=();
	$output=();

}
$i=0;

}

####################

sub extract_sff_files {

#declare var
my $name;
my $i=0;
my $output;

#declare arrays
my @filenames;
my @name;
my @name_site;

@filenames = qx(ls | grep sff);

while ($filenames[$i]) {
	$name = $filenames[$i];
	chomp $name;

	if ($name !~ /^\.+/) {
		$output = qx(mothur "#sffinfo(sff=$name,flow=F)");
	}
	$i++;
	$name=();

}
$i=0;

}

########################

sub remove_new_line_from_fasta {
	#var
	my $line;
	my $infile=${$_[0]};
	my $outfile = ${$_[1]};

	open (IN,"<",$infile) || die ("Cannot open fasta file: $!\n");
	open (TEMP,">>",$outfile) || die ("Cannot open newline removed fasta file:$!\n");

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
	my $infile = ${$_[0]};
	my $outfile = ${$_[1]};

	open (IN,"<",$infile) || die ("Can't open fasta file with newlineremoved: $!\n");
	open (TEMP,">>",$outfile) || die ("Can't open fasta file for emptyremoved$!\n");

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
	my $infile = ${$_[0]};
	my $outfile = ${$_[1]};

	#array
	my @seq;
	my @split;
	my @length;

	open (IN, "<",$infile) || die ("Error cannot read from temp1 file: $!\n");
	
	open (OUT,">>",$outfile) || die ("Error cannot write to raw stats file: $!\n");

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
	my $infile = ${$_[0]};
	my $outfile = ${$_[1]};
	
	open (IN,"<",$infile) || die ("Error opening fasta file: $!\n");
	open (TEMP,">>",$outfile) || die ("Error opening agrep formatted fasta file: $!\n");
	while (<IN>) {
		$line = $_;
		chomp $line;
		if ($line =~ /^>/) {
			$line =~ /^>(\w{14,16})\s+/; #####Edit this to work with old or new 454 IDS#####
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
	my $infile = ${$_[0]};
	my $outfile = ${$_[1]};
	my $length;
	my $line_primerremoved;

	#array
	my @line;
	my @primername;
	my @primersequence;
	my @output;
	my @fname;
	my @output2;

	open (IN,"<",$primer_file) || die "Error: $!\n";

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
		@output = qx(agrep -i -2 -D1 -I1 "^$pattern" $infile);  ##### allow 2 mismatches for a primer #####
		$fname = $primername[$i];
		#make a new folder based on first five letters of primer name
		@fname = split(//,$fname);
		$first_letter = $fname[0];
		$second_letter = $fname[1];
		$third_letter = $fname[2];
		$fourth_letter = $fname[3];
		$fifth_letter = $fname[4];
		$folder_name = $first_letter.$second_letter.$third_letter.$fourth_letter.$fifth_letter;
		# ensure first 5 letters of primer name indicate the folder results will be pooled into
		unless (-d $folder_name) {
			mkdir $folder_name;
			push(@folder_names,$folder_name);
		}

		$folder_name_path = $folder_name."/".$outfile;

		open (OUT,">>",$folder_name_path) || die "Error: $!\n";
		while ($output[$k]) {
			$line = $output[$k];
			print OUT "$line";
			$k++;
		}
		close OUT;
		$k=0;
		$i++;
		$length=();
	}

}

##########################

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
	my $dir;

	#array
	my @folders;
	my @split;
	my @output;
	my @folder_names2;

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
