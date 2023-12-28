#!/usr/bin/perl
#Oct.5, 2010 by Terri Porter
#Script to produce tabular output from clone.html files instead of the index.html file, concentrating on species level hits
#usage $ perl parse_clone_html_files.plx
#Dec.16,2010 modify to ask user to enter full path to directory containing clone files, adjust so that two summary files are provided;
#change pattern match for clone_summary_best_species.txt

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

#my $dir = "/home/terri/SAP_data/project-14.47-10.04.10_Kurtzman_fillinall_minidentity/html/clones/";
opendir (DH, $dir) || die "Error:$!\n";
@files = readdir (DH);
#print "@files\n";#test
open (OUT,">>","clone_summary_species.txt") || die ("Error cannot write to outfile:$!\n");

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
#				print "first: $id\t";#test
				print OUT "$id\t";
				$flag=1;
				}
			}
			elsif ($flag==1) {
			       if ($line =~ /\(species\)/) {
					$flag=2;
					$line =~/\+(\w+)\s{1}(\S+)\s{1}\(species\)\s{1}(\S+)%/;
					$genus = $1;
					$species = $2;
					$pp = $3;
					if ($pp =~ /~/) {
						$pp =~ s/~//;
					}
#					print "$genus\t$species\t$pp\n";#test
					print OUT "$genus\t$species\t$pp\n";
				}
				elsif ($line =~ /<p>Reliable/) {
#					print "nil\tnil\tnil\n";#test
					print OUT "nil\tnil\tnil\n";
					$flag=0;
				}
			}
			elsif ($flag==2) { 
				if ($line =~ /\(species\)/) {
					if ($line =~ /<\/a>/) {
						$line =~ />(\w+)\s{1}(\S+)<\/a>\s{1}\(species\)\s{1}(\S+)%/;
						$genus = $1;
						$species = $2;
						$pp = $3;
						if ($pp =~ /~/) {
							$pp =~ s/~//;
						}
#						print "ref: $id\t$genus\t$species\t$pp\n";#test
						print OUT "$id\t$genus\t$species\t$pp\t\n";
					}
					else {
					#print "$id\t";
					$line =~ /\+(\w+)\s{1}(\S+)\s{1}\(species\)\s{1}(\S+)%/;
					$genus = $1;
					$species = $2;
					$pp = $3;
					if ($pp =~ /~/) {
						$pp =~ s/~//;
					}
#					print "next: $id\t$genus\t$species\t$pp\n";#test
					print OUT "$id\t$genus\t$species\t$pp\n";
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

open (IN,"<","clone_summary_species.txt") || die ("Error cannot read from clone summary txt: $!\n");
@in = <IN>;
close IN;

open (OUT2,">>","clone_summary_best_species.txt") || die ("Error cannot write to clone summary best txt: $!\n");

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
	#print "line: $line\n";#test
	$line_hash{$k}=$line;
	if ($line =~ /nil\tnil\tnil\n/) {
		print OUT2 "$line\n";
		$k++;
		next;
	}
	else {
		#print "line: $line\n";#test
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
		$pp= $line_part[3];		
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


#remove species lines that contain a sp. or aff. or cf., print a tally at the end.

#var
my $z=0;
my $num_species;
my $bootstrap;
my $species2;

#array
my @infile;
my @keep;

open (IN,"<","clone_summary_best_species.txt") || die ("Error cannot read best species file: $!\n");
@infile = <IN>;
close IN;

while ($infile[$z]) {
	$line = $infile[$z];
	chomp $line;
	@line = split(/\t/,$line);
	$species2 = $line[2];
	$bootstrap = $line[3];
	if ($species2 =~ /(sp\.|cf\.|aff\.)/) {
		$z++;
		next;
	}
	else {
		if ($bootstrap >= 95) {
			push(@keep,$line);
		}
	}
	$z++;
}
$num_species = scalar(@keep);
print "There are $num_species\n";
