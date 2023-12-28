#!/usr/bin/perl
#Oct.5, 2010 by Terri Porter
#Script to produce tabular output from clone.html files instead of the index.html file, concentrating on species level hits
#usage $ perl parse_clone_html_files.plx
#Dec.16,2010 modify to ask user to enter full path to directory containing clone files, adjust so that two summary files are provided;
#edit pattern to grab appropriate gi number from query seq fasta style header

use strict;
use warnings;

#declare var
my $dir;
my $pattern;
my $file;
my $line;
my $id;
my $i=0;
my $j=0;
my $flag=0;
my $genus;
my $species;
my $pp;

#declare array
my @files;
my @current_file;
my @line;

print "Enter full path to directory containing clone files (must include trailing '/'):\n";
$dir = <STDIN>;
chomp $dir;

opendir (DH, $dir) || die "Error:$!\n";
@files = readdir (DH);

open (OUT,">>","clone_summary_genera.txt") || die ("Error cannot write to outfile:$!\n");

open (TRASH,">>","trash.txt") || die ("Error cannot open trash file :$!\n");

print "Enter basename for clone files:\n";
$pattern = <STDIN>;
chomp $pattern;

while ($files[$i]) {
	$file = $files[$i];
	if ($file =~ /$pattern/) {

		open (FH, $file) || die ("Error:$!\n");
		@current_file = <FH>;
		while ($current_file[$j]){
			$line = $current_file[$j];
			chomp $line;
			if ($flag==0) {
		       		if ($line =~ /<h1>Sequence:/) {
				$line =~ /<h1>Sequence: (.+)<\/h1>/g;
				$id = $1;

				print OUT "$id\t";
				$flag=1;
				}
				#		else {
				#	print TRASH "$file\n";#test find entries with no sequence
				#}
			}
			elsif ($flag==1) {
			       if ($line =~ /\(genus\)/) {

				       if ($line =~/<\/a>/) {
						$flag=2;
						$line =~/>(\w+)<\/a>\s{1}\(genus\)\s{1}(\S+)%/;
						$genus = $1;
						$pp = $2;
						if ($pp =~ /~/) {
							$pp =~ s/~//;
						}
						print OUT "$genus\t$pp\n";
					}
					else {
						$flag=2;
						$line =~ /\+(\w+)\s{1}\(genus\)\s{1}(\S+)%/;
						$genus = $1;
						$pp = $2;
						if ($pp =~ /~/) {
							$pp =~ s/~//;
						}
						print OUT "$genus\t$pp\n";
					}
				}
				elsif ($line =~ /<p>Reliable/) {
					print OUT "nil\tnil\n";
					print TRASH "$file\n";#get entries with no classification
					$flag=0;
				}
			}
			elsif ($flag==2) { 
				if ($line =~ /\(genus\)/) {
					if ($line =~ /<\/a>/) {
						$line =~ />(\w+)<\/a>\s{1}\(genus\)\s{1}(\S+)%/;
						$genus = $1;
						$pp = $2;
						if ($pp =~ /~/) {
							$pp =~ s/~//;
						}
						print OUT "$id\t$genus\t$pp\t\n";
					}
					else {
						$line =~ /\+(\w+)\s{1}\(genus\)\s{1}(\S+)%/;
						$genus = $1;
						$pp = $2;
						if ($pp =~ /~/) {
							$pp =~ s/~//;
						}
						print OUT "$id\t$genus\t$pp\n";
					}
				}
				elsif ($line =~ /<\/pre><\/div>/) {
					$flag=0;
				}
				else {
					$j++;
					next;
				}
			}
			$j++;
		}
		$j=0;
		$i++;
		close FH;
	}
	else {
		$i++;
		next;
	}
}
close OUT;

#array
my @in;
my @values;
my @keys;
my @pp_part;
@line=();
my @line_part;

open (IN,"<","clone_summary_genera.txt") || die ("Error cannot read from clone summary txt: $!\n");
@in = <IN>;
close IN;

open (OUT2,">>","clone_summary_best_genera.txt") || die ("Error cannot write to clone summary best txt: $!\n");

#var
my $k=0;
my $current_line;
my $id_part;
my $id_current;
my $id_previous="nil";
my $flag2 = 0;
my $key;
my $pp_value;
my $best_pp;
my $line_index;
my $best_line;
my $pp_part;

#hash
my %pp_hash;
my %line_hash;

while ($in[$k]) {
	$line = $in[$k];
	chomp $line;
	print "line: $line\n";#test
	$line_hash{$k}=$line;
	if ($line =~ /nil\tnil\n/g) {
		print OUT2 "$line\n";
		$k++;
		next;
	}
	else {
		print "line2: $line\n";#test
		$current_line = $line_hash{$k};
		print "current_line: $current_line\n";
		@line_part = split(/\t/,$current_line);
		$id_part = $line_part[0];
	        print "id_part: $id_part\n";#test
		#$id_part =~ /\d+\|\*\|(\w{14})/;###edit pattern here!###
		my @id_part = split(/\|/,$id_part);
       		$id_current = $id_part[1];
		print "id_current: $id_current\n";#test
		#print "arrayline: @line\n";#test
		$pp= pop(@line_part);		
		#print "pp part: $pp_part\n";#test
		#@pp_part=split(/\t/, $pp_part);
		#$pp = pop(@pp_part);
		print "pp: $pp\n";#test
		if ($id_current ne $id_previous) {
			if ($flag2==0)	{
				$pp_hash{$k}=$pp;
				$id_previous = $id_current;
				$flag2=1;
			}
			elsif ($flag2==1) {
				
				foreach $key (sort{$pp_hash{$b}<=>$pp_hash{$a}} keys %pp_hash) {
					$pp_value=$pp_hash{$key};
					push(@values,$pp_value);
					push(@keys,$key);
				}
				$best_pp = $values[0];
				$line_index = $keys[0];
				$best_line = $line_hash{$line_index};
				print OUT2 "$best_line\n";
				@values=();
				@keys=();
				%pp_hash=();

				$pp_hash{$k}=$pp;
				$id_previous = $id_current;
			}
		}
		elsif ($id_current eq $id_previous) {
			$pp_hash{$k}=$pp;
			$id_previous = $id_current;
		}
	}
	$k++;
}

#don't forget to parse last set of entries!

foreach $key (sort{$pp_hash{$b}<=>$pp_hash{$a}} keys %pp_hash) {
	$pp_value=$pp_hash{$key};
	push(@values,$pp_value);
	push(@keys,$key);
}
$best_pp = $values[0];
$line_index = $keys[0];
$best_line = $line_hash{$line_index};
print OUT2 "$best_line\n";
close OUT2;

#keep only best results with >80% ConstrainedNJ bootstrap proportion

#var
my $bootstrap;
my $z=0;
my $num_keep;

#array
my @infile;
my @line;
my @keep;

open (IN,"<","clone_summary_best_genera.txt") || die ("Error can't read best genera file: $!\n");
@infile = <IN>;
close IN;

while ($infile[$z]) {
	$line = $infile[$z];
	chomp $line;
	@line = split(/\t/,$line);
	$bootstrap = $line[2];

	if ($bootstrap >= 95) {
		push(@keep,$line);
	}
	else {
		$z++;
		next;
	}
	$z++;
}

$num_keep = scalar(@keep);

print "There are $num_keep genera with support >= 95\n";
	
