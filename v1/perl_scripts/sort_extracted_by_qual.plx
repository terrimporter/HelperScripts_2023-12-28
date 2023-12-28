#!/usr/bin/perl
#July 8, 2011 by Terri Porter
#Script to sort extracted ITS sequences by average phred quality before running usearch
#usage perl sort_extracted_by_qual.plx outdata.csv seqtrim.fasta.newlineremoved seqtrim.qual.newlineremoved
 
use strict;
use warnings;

#declare var
my $key;
my $value;
my $i=0;
my $line;

#declare array
my @in;
my @seqtrim_fasta;
my @extracted;
my @ITS1;
my @ITS2;
my @reverse_complemented;
my @ITS1_filtered;
my @ITS2_filtered;
my @qual;
my @phred_extracted;
my @phred_filtered;
my @ITS1_phred_filtered;
my @just_reversed;
my @ITS2_phred_filtered;
my @fasta;
my @phred;
my @sorted_fasta;
my @sorted_qual;
my @sorted_fasta_in;
my @sorted_qual_in;
my @qual_sorted_out;
my @fasta_sorted_out;

#declare hash
my %ITS1_start;
my %ITS1_end;
my %ITS2_start;
my %ITS2_end;
my %ITS_revcomp;
my %ITS1_seq;
my %ITS2_seq;
my %ITS1_phred_seq;
my %ITS2_phred_seq;
my %fasta;

open (IN,"<",$ARGV[0]) || 
	die ("Error reading infile: $!\n");
	@in = <IN>;	
close IN;

open (IN2,"<",$ARGV[1]) || 
	die ("Error reading seqtrim.fasta: $!\n");
	@seqtrim_fasta = <IN2>;
close IN2;

get_extraction_range();

extract_fasta();

open (TMP1, ">>", "ITS1.extracted.temp") || 
	die ("Error writing ITS1.temp: $!\n");
open (TMP2, ">>", "ITS2.extracted.temp") || 
	die ("Error writing ITS2.temp: $!\n");

	while ( ($key,$value) = each (%ITS1_seq)) {
		print TMP1 ">$key\n$value\n";
	}

	while ( ($key,$value) = each (%ITS2_seq)){
		print TMP2 ">$key\n$value\n";
	}
close TMP1;
close TMP2;

open (IN, "<","ITS1.extracted.temp") || 
	die ("Error reading ITS1.temp: $!\n");
	@ITS1 = <IN>;
close IN;

filter_fasta_by_length(@ITS1);

open (TMP3,">>","ITS1.extracted.filtered.temp") || 
	die ("Error writing ITS1.extracted.filtered.temp: $!\n");

	while ($extracted[$i]) {
		$line = $extracted[$i];
		print TMP3 "$line\n";
		$i++;
	}
	$i=0;
close TMP3;

open (IN2, "<","ITS2.extracted.temp") || 
	die ("Error reading ITS2.temp: $!|n");
	@ITS2 = <IN2>;
close IN2;

filter_fasta_by_length(@ITS2);

open (TMP4,">>","ITS2.extracted.filtered.temp") || 
	die ("Error writing ITS2.extracted.filtered.temp: $!\n");

	while ($extracted[$i]) {
		$line = $extracted[$i];
		print TMP4 "$line\n";
		$i++;
	}
	$i=0;
close TMP4;

open (IN3,"<","ITS1.extracted.filtered.temp") || 
	die ("Error reading ITS1.extracted.filtered.temp: $!\n");
	@ITS1_filtered = <IN3>;
close IN3;

rev_comp(@ITS1_filtered);

open (TMP5,">>","ITS1.extracted.filtered.rc.temp") || 
	die ("Error reading ITS1.extracted.filtered.rc.temp: $!\n");

	while ($reverse_complemented[$i]) {
		$line = $reverse_complemented[$i];
		print TMP5 "$line\n";
		$i++;
	}
	$i=0;
close TMP5;

open (IN4,"<","ITS2.extracted.filtered.temp") || 
	die ("Error reading ITS2.extracted.filtered.temp: $!\n");
	@ITS2_filtered = <IN4>;
close IN4;

rev_comp(@ITS2_filtered);

open (TMP6,">>","ITS2.extracted.filtered.rc.temp") || 
	die ("Error reading ITS2.extracted.filtered.rc.temp: $!\n");

	while ($reverse_complemented[$i]) {
		$line = $reverse_complemented[$i];
		print TMP6 "$line\n";
		$i++;
	}
	$i=0;
close TMP6;

open (IN5,"<",$ARGV[2]) || 
	die ("Error cannot read infile.qual: $!\n");
	@qual = <IN5>;
close IN5;

extract_phred(@qual);

open (TMP7,">>","ITS1.phred.extracted.temp") || 
	die ("Error writing ITS1.phred.extracted.temp:$!\n");

	while ( ($key,$value) = each (%ITS1_phred_seq)){
		print TMP7 ">$key\n$value\n";
	}	       
close TMP7;

open (TMP8, ">>","ITS2.phred.extracted.temp") || 
	die ("Error writing ITS2.phred.extracted.temp: $!\n");

	while ( ($key,$value) = each (%ITS2_phred_seq)) {
		print TMP8 ">$key\n$value\n";
	}
close TMP8;

open (IN6,"<","ITS1.phred.extracted.temp") || 
	die ("Error reading ITS1.phred.extracted.temp:$!\n");
	@phred_extracted=<IN6>;
close IN6;

filter_phred_by_length(@phred_extracted);

open (TMP9,">>","ITS1.phred.extracted.filtered.temp") || 
	die ("Error writing ITS1.phred.extracted.filtered.temp: $!\n");

	while ($phred_filtered[$i]) {
		$line = $phred_filtered[$i];
		print TMP9 "$line\n";
		$i++;
	}
	$i=0;
close TMP9;

open (IN7,"<","ITS2.phred.extracted.temp") || 
	die ("Error reading ITS2.phred.extracted.temp: $!\n");
	@phred_extracted=<IN7>;
close IN7;

filter_phred_by_length(@phred_extracted);

open (TMP10,">>","ITS2.phred.extracted.filtered.temp") || 
	die ("Error writing ITS2.phred.extracted.filtered.temp:$!\n");

	while ($phred_filtered[$i]) {
		$line = $phred_filtered[$i];
		print TMP10 "$line\n";
		$i++;
	}
	$i=0;
close TMP10;

open (IN8,"<","ITS1.phred.extracted.filtered.temp") || 
	die ("Error reading ITS1.phred.extracted.filtered.temp:$!\n");
	@ITS1_phred_filtered = <IN8>;
close IN8;

rev(@ITS1_phred_filtered);

open (TMP11,">>","ITS1.phred.extracted.filtered.rev.temp") || 
	die ("Error writing ITS1.phred.extracted.filtered.rev.temp:$!\n");

	while ($just_reversed[$i]) {
		$line = $just_reversed[$i];
		print TMP11 "$line\n";
		$i++;
	}
	$i=0;
close TMP11;

open (IN9,"<","ITS2.phred.extracted.filtered.temp") ||
	die ("Error reading ITS2.phred.extracted.filtered.temp: $!\n");
	@ITS2_phred_filtered = <IN9>;
close IN9;

rev(@ITS2_phred_filtered);

open (TMP12,">>","ITS2.phred.extracted.filtered.rev.temp") ||
	die ("Error writing ITS2.phred.extracted.filtered.rev.temp:$!\n");

	while ($just_reversed[$i]) {
		$line = $just_reversed[$i];
		print TMP12 "$line\n";
		$i++;
	}
	$i=0;
close TMP12;

open (IN10,"<","ITS1.extracted.filtered.rc.temp") ||
	die ("Error reading ITS1.extracted.filtered.rc.temp:$!\n");
	@fasta = <IN10>;
close IN10;

open (IN11,"<","ITS1.phred.extracted.filtered.rev.temp") ||
	die ("Error reading ITS1.phred.extracted.filtered.rev.temp:$!\n");
	@phred = <IN11>;
close IN11;

sort_by_qual(\@fasta,\@phred); #pass two array refs to subroutine

open (TMP13,">>","ITS1_sorted_by_qual.fasta.temp") ||
	die ("Error cannot write to ITS1_sorted_by_qual.temp:$!\n");

	while ($sorted_fasta[$i]) {
		$line = $sorted_fasta[$i];
		print TMP13 "$line\n";
		$i++;
	}
	$i=0;
close TMP13;

open (TMP14,">>","ITS1_sorted_by_qual.qual.temp") ||
	die ("Error cannot write to ITS1_sorted_by_qual.temp: $!\n");

	while ($sorted_qual[$i]) {
		$line = $sorted_qual[$i];
		print TMP14 "$line\n";
		$i++;
	}
	$i=0;
close TMP14;

open (IN12,"<","ITS2.extracted.filtered.rc.temp") ||
	die ("Error reading ITS2.extracted.filtered.rc.temp:$!\n");
	@fasta = <IN12>;
close IN12;

open (IN13,"<","ITS2.phred.extracted.filtered.rev.temp") ||
	die ("Error reading ITS2.phred.extracted.filtered.rev.temp:$!\n");
	@phred = <IN13>;
close IN13;

sort_by_qual(\@fasta,\@phred);

open (TMP15,">>","ITS2_sorted_by_qual.fasta.temp") ||
	die ("Error cannot write to ITS2_sorted_by_qual.temp:$!\n");

	while ($sorted_fasta[$i]) {
		$line = $sorted_fasta[$i];
		print TMP15 "$line\n";
		$i++;
	}
	$i=0;
close TMP15;

open (TMP16,">>","ITS2_sorted_by_qual.qual.temp") ||
	die ("Error cannot write to ITS2_sorted_by_qual.temp:$!\n");

	while ($sorted_qual[$i]) {
		$line = $sorted_qual[$i];
		print TMP16 "$line\n";
		$i++;
	}
	$i=0;
close TMP16;

open (IN14,"<","ITS1_sorted_by_qual.fasta.temp") ||
	die ("Error cannot read ITS1_sorted_by_qual.fasta.temp:$!\n");
	@sorted_fasta_in = <IN14>;
close IN14;

open (IN15,"<","ITS1_sorted_by_qual.qual.temp") ||
	die ("Error cannot read ITS1_sorted_by_qual.qual.temp:$!\n");
	@sorted_qual_in = <IN15>;
close IN15;

sort_by_qual_for_usearch(\@sorted_fasta_in, \@sorted_qual_in);

open (OUT1,">>","ITS1_extracted_sorted.fasta") ||
	die ("Error cannot write to ITS1_extracted_sorted.fasta:$!\n");

	while ($fasta_sorted_out[$i]) {
	       $line = $fasta_sorted_out[$i];
		chomp $line;
		print OUT1 "$line\n";
		 $i++;
	}
	$i=0;
close OUT1;

open (OUT2,">>","ITS1_extracted_sorted.qual") ||
	die ("Error cannot write to ITS1_extracted_sorted.qual:$!\n");

	while ($qual_sorted_out[$i]) {
		$line = $qual_sorted_out[$i];
		chomp $line;
		print OUT2 "$line\n";
		$i++;
	}
	$i=0;
close OUT2;

open (IN16,"<","ITS2_sorted_by_qual.fasta.temp") ||
	die ("Error cannot read ITS2_sorted_by_qual.fasta.temp:$!\n");
	@sorted_fasta_in = <IN16>;
close IN16;

open (IN17,"<","ITS2_sorted_by_qual.qual.temp") ||
	die ("Error cannot read ITS2_sorted_by_qual.qual.temp: $!\n");
	@sorted_qual_in = <IN17>;
close IN17;

sort_by_qual_for_usearch(\@sorted_fasta_in, \@sorted_qual_in);

open (OUT3,">>","ITS2_extracted_sorted.fasta") ||
	die ("Error cannot write to ITS2_extracted_sorted.fasta:$!\n");

	while ($fasta_sorted_out[$i]) {
		$line = $fasta_sorted_out[$i];
		chomp $line;
		print OUT3 "$line\n";
		$i++;
	}
	$i=0;
close OUT3;

open (OUT4,">>","ITS2_extracted_sorted.qual") ||
	die ("Error cannot write to ITS2_extracted_sorted.qua:$!\n");

	while ($qual_sorted_out[$i]) {
		$line = $qual_sorted_out[$i];
		chomp $line;
		print OUT4 "$line\n";
		$i++;
	}
	$i=0;
close OUT4;

####################

sub get_extraction_range {

#declare var
my $i=0;
my $line;
my $id;
my $start;
my $end;
my $start2;
my $end2;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /\)\s+\w{14}\w+/) {
		$line =~ /\)\s+(\w{14})\w+/;
		$id = $1;
		
		if ($line =~ /ITS1:\s{1}\w+-\w+\s+/) {
			$line =~ /ITS1:\s{1}(\w+)-(\w+)\s+/;
			$start = $1;
			$end = $2;
			$ITS1_start{$id} = $start;
			$ITS1_end{$id} = $end;
		}
		if ($line =~ /ITS2:\s{1}\w+-\w+/) {
			$line =~ /ITS2:\s{1}(\w+)-(\w+)/;
			$start2 = $1;
			$end2 = $2;
			$ITS2_start{$id} = $start2;
			$ITS2_end{$id} = $end2;
		}
		if ($line =~ /Reverse complementary/) {
			$ITS_revcomp{$id}=1;
		}
	}
	$i++;
}

}

####################

sub extract_fasta {

#declare var
my $i=0;
my $line;
my $seqtrim_id;
my $start;
my $end;
my $j;
my $next_line;
my $length;
my $start_substring;
my $length_substring;
my $ITS1_seq;
my $start2;
my $end2;
my $next_line2;
my $length2;
my $start_substring2;
my $length_substring2;
my $ITS2_seq;
my $k=0;

#declare array
my @sequence;
my @sequence2;

while ($seqtrim_fasta[$i]) {
	$line =$seqtrim_fasta[$i];
	chomp $line;

	if ($line =~ /^>\w{14}\s+/) {
		$line =~ /^>(\w{14})\s+/;
		$seqtrim_id = $1;

		if ($ITS1_start{$seqtrim_id}) {
			$start = $ITS1_start{$seqtrim_id};
			$end = $ITS1_end{$seqtrim_id};

			$j=$i+1;
			$next_line = $seqtrim_fasta[$j];
			chomp $next_line;
			@sequence = split (//,$next_line);
			$length = scalar(@sequence);

			if ($end eq 'end') {
				$end = $length;
			}
			
			$start_substring = $start-1;
			$length_substring = $end-$start;
			$ITS1_seq = substr($next_line, $start_substring, $length_substring);
			$ITS1_seq{$seqtrim_id} = $ITS1_seq;
					
			@sequence=();
			$ITS1_seq=();
		}

		if ($ITS2_start{$seqtrim_id}) {
			$start2 = $ITS2_start{$seqtrim_id};
			$end2 = $ITS2_end{$seqtrim_id};

			$j=$i+1;
			$next_line2 = $seqtrim_fasta[$j];
			chomp $next_line2;

			@sequence2 = split(//,$next_line);
			$length2 = scalar(@sequence2);

			if ($end2 eq 'end') {
				$end2 = $length2;
			}

			$start_substring2 = $start2-1;
			$length_substring2 = $end2-$start2;
			$ITS2_seq = substr($next_line2, $start_substring2, $length_substring2);
			$ITS2_seq{$seqtrim_id} = $ITS2_seq;

			@sequence2=();
			$ITS2_seq=();
		}
	}
	$i++;
}

}

####################

sub filter_fasta_by_length {

#declare var
my $i=0;
my $line;
my $header;
my $j;
my $next_line;
my $length;
my $sequence;

#declare array
my @next_line;
my @fasta=();

@fasta = @_;
#print "@fasta";

@extracted=();

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		#print $header."\n";#test
		$j=$i+1;
		$next_line = $fasta[$j];
		chomp $next_line;
		@next_line = split(//,$next_line);
		$length = scalar(@next_line);

		if ($length >= 100) {##### change minimum sequence length here #####
			$sequence = join('',@next_line);
			push(@extracted,$header);
			push(@extracted,$sequence);
			$sequence=();
		}
		$i+=2;
		@next_line=();
	}
	else {
		$i++;
	}
}

}

####################

sub rev_comp {

#declare var
my $i=0;
my $line;
my $id;
my $j;
my $next_line;
my $k=0;
my $base;
my $reverse_complemented;

#declare array
my @next_line;
my @reversed;
my @comp=();
my @fasta=();
@reverse_complemented=();

@fasta = @_;

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ /^>(\w{14})/;
		$id = $1;

		if ($ITS_revcomp{$id}) {
			$j=$i+1;
			$next_line = $fasta[$j];
			chomp $next_line;
			@next_line = split(//,$next_line);
			@reversed = reverse(@next_line); #reverse sequence

			while ($reversed[$k]) {
				$base = $reversed[$k];
				#complement bases
				if (($base eq 'A') || ($base eq 'a')) {
					push(@comp, 'T');
				}
				elsif (($base eq 'C') || ($base eq 'c')) {
					push (@comp, 'G');
				}
				elsif (($base eq 'T') || ($base eq 't')) {
					push (@comp,'A');
				}
				elsif (($base eq 'G') || ($base eq 'g')) {
					push (@comp, 'C');
				}
				elsif (($base eq 'N') || ($base eq 'n')) {
					push (@comp ,'N');
				}
				$k++;
			}
			$k=0;
			$reverse_complemented = join("",@comp);
			push(@reverse_complemented, $line);
			push(@reverse_complemented, $reverse_complemented);
			@comp=();
		}
		else {
			push(@reverse_complemented, $line);
			$j=$i+1;
			$next_line = $fasta[$j];
			chomp $next_line;
			push(@reverse_complemented, $next_line);
		}
		$i+=2;
	}
	else {
		$i++;
	}
}
	
}

####################

sub extract_phred {

#declare var
my $i=0;
my $line;
my $seqtrim_id;
my $start;
my $end;
my $j;
my $phred_line;
my $length;
my $start_substring;
my $end_substring;
my $k=0;
my $score;
my $ITS1_phred_seq=();
my $start2;
my $end2;
my $start_substring2;
my $end_substring2;
my $ITS2_phred_seq=();

#declare array
my @seqtrim_qual=();
my @phred_sequence=();
my @phred_subsequence;

@seqtrim_qual = @_;

while ($seqtrim_qual[$i]) {
	$line =$seqtrim_qual[$i];
	chomp $line;

	if ($line =~ /^>\w{14}\s+/) {
		$line =~ /^>(\w{14})\s+/;
		$seqtrim_id = $1;

		if ($ITS1_start{$seqtrim_id}) {
			$start = $ITS1_start{$seqtrim_id};
			$end = $ITS1_end{$seqtrim_id};
		
			$j=$i+1;
                        $phred_line = $seqtrim_qual[$j];
                        chomp $phred_line;
                        @phred_sequence = split (/ /,$phred_line);
			$length = scalar(@phred_sequence);
			
			if ($end eq 'end') {
                        	$end = $length;
                        }

			$start_substring = $start-1;
			$end_substring = $end-1;
			while ($phred_sequence[$k]) {
				$score = $phred_sequence[$k];
				if ($k>= $start_substring) {
					if ($k<=$end_substring) {
						push(@phred_subsequence,$score);
					}
				}
				$k++;
			}
			$k=0;
			$ITS1_phred_seq = join(" ",@phred_subsequence);
			$ITS1_phred_seq{$seqtrim_id} = $ITS1_phred_seq;
												 
			 @phred_sequence=();
			 @phred_subsequence=();
                         $ITS1_phred_seq=();
		}

		if ($ITS2_start{$seqtrim_id}) {
			$start2 = $ITS2_start{$seqtrim_id};
			$end2 = $ITS2_end{$seqtrim_id};

			$j=$i+1;
			$phred_line = $seqtrim_qual[$j];
			chomp $phred_line;
			@phred_sequence = split (/ /,$phred_line);
			$length = scalar(@phred_sequence);
			
			if ($end2 eq 'end') {
                        	$end2 = $length;
			}

			$start_substring2 = $start2-1;
			$end_substring2 = $end2-1;
			while ($phred_sequence[$k]) {
				$score = $phred_sequence[$k];
				if ($k>= $start_substring2) {
					if ($k<=$end_substring2) {
						push(@phred_subsequence,$score);
					}
				}
				$k++;
			}
			$k=0;
		
			$ITS2_phred_seq = join(" ",@phred_subsequence);
		        $ITS2_phred_seq{$seqtrim_id} = $ITS2_phred_seq;	

			@phred_sequence=();
			@phred_subsequence=();
			$ITS2_phred_seq=();
		}		
		$i+=2;
	}
	else {
      		$i++;
	}				
}

}

####################

sub filter_phred_by_length {

#declare var
my $i=0;
my $line;
my $header;
my $id;
my $j;
my $length;
my $phred_line;

#declare array
my @phred_line;
my @qual=();

@phred_filtered=();
@qual = @_;

while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w{14})/;
		$id = $1;

		$j=$i+1;
		$phred_line = $qual[$j];
		chomp $phred_line;

		@phred_line = split(/ /,$phred_line);
		$length = scalar(@phred_line);

		if ($length >= 100) {##### change minimum sequence length here #####
			push(@phred_filtered,$header);
			push(@phred_filtered,$phred_line);
			$phred_line=();
			@phred_line=();
		}
		$i+=2;
	}
	else {
		$i++;
	}
}

}

####################

sub rev {

#declare var
my $line;
my $id;
my $i=0;
my $j;
my $phred_line;
my $reversed;

#declare array
my @phred_fasta=();
my @phred_line;
my @reversed;

@just_reversed=();
@phred_fasta = @_;

while ($phred_fasta[$i]) {
	$line = $phred_fasta[$i];
	chomp $line;
	
	if ($line =~ /^>/) {
		$line =~ /^>(\w{14})/;
		$id = $1;
		#print $id."\n";#test
		if ($ITS_revcomp{$id}) {
			$j=$i+1;
			$phred_line = $phred_fasta[$j];
			chomp $phred_line;
			@phred_line = split(/ /,$phred_line);
			@reversed = reverse(@phred_line); #reverse sequence
			$reversed = join(" ",@reversed);
			#print "found rev comp in hash\n";#test
			push(@just_reversed, $line);
			push(@just_reversed, $reversed);
		}
		else {
			#print "didn't find id in revcomp hash\n";#test
			push(@just_reversed, $line);
			$j=$i+1;
			$phred_line = $phred_fasta[$j];
			chomp $phred_line;
			push(@just_reversed, $phred_line);
		}
		$i+=2;
	}
	else {
		$i++;
	}
}

}

####################

sub sort_by_qual {

#declare var
my $line;
my $header;
my $id;
my $i=0;
my $j=0;
my $sequence;
my $k;
my $phred_header;
my $phred_id;
my $phred_sequence;

#declare array
my @fasta;
my @phred;
my @ids;

@sorted_fasta=();
@sorted_qual=();
@fasta = @{$_[0]}; #dereference
@phred = @{$_[1]};

while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w{14})/;
		$id = $1;
		$j=$i+1;
		$sequence = $fasta[$j];
		chomp $sequence;
		$fasta{$id} = $sequence;
		$i+=2;
	}
	else {
		$i++;
	}
}
$i=0;

while ($phred[$i]) {
	$line = $phred[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$phred_header = $line;
	 	$phred_header =~ /^>(\w{14})/;
		$phred_id = $1;

		if ($fasta{$phred_id}) {
			$header = ">".$phred_id;
			$sequence = $fasta{$phred_id};
			push(@sorted_fasta,$header);
			push(@sorted_fasta,$sequence);
			push(@sorted_qual,$phred_header);
			$j = $i+1;
			$phred_sequence = $phred[$j];
			chomp $phred_sequence;
			#print "$phred_sequence\n";#test
			push(@sorted_qual,$phred_sequence);
			
		}
		$i+=2;
	}
	else {
		$i++;
	}
}

}

####################

sub sort_by_qual_for_usearch {

#declare var
my $i=0;
my $line;
my $j=0;
my $phred_score;
my $k=0;
my $l;
my $index;
my $header;
my $phred_seq;
my $m=0;
my $fasta_header;
my $fasta_seq;

#declare array
my @fasta;
my @qual;
my @phred_header=();
my @phred_seq=();
my @phred_scores=();
my @value;
my @fasta_header;
my @fasta_seq;

#declare hash
my %hash=();

#dereference
@fasta = @{$_[0]};
@qual = @{$_[1]};


@qual_sorted_out=();
@fasta_sorted_out=();

while($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		push (@phred_header, $line);
	}
	else {
		push (@phred_seq, $line);
	}
	$i++;
}

while ($phred_seq[$j]) {
	$line = $phred_seq[$j];
	@phred_scores = split(/ /,$line);
	foreach $phred_score (@phred_scores) {
		if ($phred_score < 20) {
			$k++;
		}
	}
	$hash{$j} = $k;
	@phred_scores=();
	$k=0;
	$j++;
}

foreach $phred_score (sort{$hash{$a} <=> $hash{$b}} keys %hash) {
	push (@value, $phred_score);
}

foreach $l (@value) {
	$index = $l;
	$header = $phred_header[$index];
	$phred_seq = $phred_seq[$index];
	push(@qual_sorted_out,$header);
	push(@qual_sorted_out,$phred_seq);
}

while($fasta[$m]) {
	$line = $fasta[$m];
	chomp $line;
	if ($line =~ /^>/) {
		push(@fasta_header,$line);
	}
	else {
		push(@fasta_seq,$line);
	}
	$m++;
}

foreach $l (@value) {
	$index = $l;
	$fasta_header = $fasta_header[$index];
	$fasta_seq = $fasta_seq[$index];
	push(@fasta_sorted_out,$fasta_header);
	push(@fasta_sorted_out,$fasta_seq);
}

}
