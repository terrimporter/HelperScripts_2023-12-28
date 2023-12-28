#!/usr/bin/perl
#August 21, 2012 edited to check for . and .. files differently using /^\/ instead of pop-ing files off array
#September 30, 2009
#Script to reformat gblocks.fasta by truncating taxon names to the first two characters for consistency among alignments, removing gaps from sequence alignments, rearranging old alignment so that it's in alphabetical order, inserting missing taxa with gapped sequences where necessary
#USAGE $perl reformat_gblocks.plx relative/path/to/gblocks/infiles absolute/path/to/formatted_gblocks/outfiles

use strict;
use warnings;

#declare some counting variables
my $dir;
my $i=0;#index file names
my $infile;
my $line;
my $truncated;
my $j=0;#flags continuation of previous sequence line
my $seq1;
my $seq2;
my $seq_length;
my $longest_seq_length;
my $flag=0;#flags files with no gblock sequence data
my $count_gap=1;
my $flag2=0;#flags element in @header not in the same order as @ref elements
my $gap_seq;
my $header;
my $value;
my $key;
my $count_old=0;#index @header and @ref elements

#declare some arrays
my @gblock_files;
my @header;
my @seq;
my @split;
my @count_seq_length;
my @count_seq_length_sorted;
my @gap;

#declare some hashes
my %aln;

#read a directory of file names into an array
$dir = "$ARGV[0]";
opendir DH, $dir;
@gblock_files = readdir(DH);

#at unix prompt, do $ls -a to see if the . and .. files are present, if so then:
#pop(@gblock_files);#to remove .
#pop(@gblock_files);#to remove ..

while ($gblock_files[$i]){
	#open infile and prepare outfile
	$infile = $ARGV[0].$gblock_files[$i]; #path/filename
	if ($infile !~ /^\./) {
		open (IN, "<", $infile) || die "Can't read $infile: $!";

		#read in fasta file and store data in two arrays, one for fasta headers, one for sequences
		while (<IN>){
			$line = $_;
			chomp($line);
			if (/>/ && $j==0){#first fasta header
				$truncated = substr $line, 1, 2;#offset by 1 to avoid >, keep first two char of header
				push (@header, $truncated);#add to end of array
			}
			unless (/>/) {#sequence line
				$line =~ s/\s+//g;#removes any whitespace
				if ($j==0) {
					$seq1 = $line;
					$j=1;
				}
				elsif ($j==1) {
				$seq2 = $line;
				$seq1 = $seq1.$seq2;#successively concatenate each line of sequence
				}
			}	
			elsif (/>/ && $j==1) {#fasta header
				push (@seq,$seq1);#add to end of array
				$j=0;#reset flag
				$truncated = substr $line,1,2;
				push (@header, $truncated);#add to end of array
			}	
		}
		$j=0;#reset flag
		push (@seq, $seq1);#be sure to add last sequence of file to array!
	
		foreach(@seq) {
			$line = $_;
			if (/\S+/) {#search for non-white space in sequence line
				#count length of sequence
				@split = split(//,$line);
				$seq_length = scalar(@split);
				push (@count_seq_length, $seq_length);#add to end of array
			}
			else {
				$seq_length  = 0;
				push (@count_seq_length, $seq_length);#add to end of array
			}
		}		
		@count_seq_length_sorted = sort {$b <=> $a} @count_seq_length;#numerical descending sort
		$longest_seq_length = shift(@count_seq_length_sorted);#remove first element off array which should be length of longest sequence
	
		if ($longest_seq_length > 0){
			$flag=0;
		}
		if ($longest_seq_length == 0) {#flags files with no gblock sequence data
			$flag=1;
		}

		#only process files with gblock sequence data
		if ($flag==0) {
			#create a new sequence of gap characters of appropriate length
			while ($count_gap <= $longest_seq_length){
				push (@gap, "-");#add to end of array
				$count_gap++;
			}
		
			$count_gap = 1;#reset counter for next file
			$gap_seq = join("",@gap);

			#define a hash with keys equivalent to all taxon ids
			%aln = ("AM"=> '',"AN"=> '',"BD"=> '',"BE"=> '',"GP"=> '',"LB"=> '',"MB"=> '',"MV"=> '',"MC"=> '',"PB"=> '',"PE"=> '',"RO"=> '',"SC"=> '',"SP"=> '',"SQ"=> '',"SR"=> '',"UM"=> '',"AF"=> '', "CR"=> '', "NE" => '', "PD"=> '', "EP"=> '', "BB"=> '', "CA"=> '');

			#for each hash key, check against @header, if match found, assign hash value from @ref
			for $key (sort keys %aln){
				foreach $header (@header){
					if ($key eq $header) {
						$aln{$key} = $seq[$count_old];#assign seq to appropriate taxid in hash
					}
					$count_old++;
				}
				$count_old=0;#reset counter
			}

			#add gapped sequences for missing taxa
			while (($key,$value)=each(%aln)){
				if ($aln{$key} !~ /\S+/) {#search for lack of characters
					$aln{$key}=$gap_seq;
				}
			}

			my $outfile = $ARGV[1].$gblock_files[$i]."\.fasta"; #generate outfile name from infile, adding .fasta extension
    	    open (OUT, ">>", $outfile) || die "Can't read $outfile: $!";
	
			#sort keys then print OUT
			foreach $key (sort keys %aln) {
				print OUT $key."\t".$aln{$key}."\n";
			}
		}
		@header=();#empty array
		@seq=();
		@split=();
		@count_seq_length=();
		@count_seq_length_sorted=();
		@gap=();
		%aln=();#empty hash
		close IN;
		close OUT;
	}
	$i++;
}
