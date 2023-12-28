#!/usr/bin/perl
#Oct.20,2011 by Terri Porter
#Script to use list of gi's from get_gi_from_fasta_header.plx
#and grab corresponding blastn files and write them to a merged file
#usage perl grab_blastn_from_dir.plx gi.txt

#declare var
my $dir;
my $file;
my $gi;
my $i=0;
my $gi_ref;
my $path_to_file;
my $line;

#declare array
my @gi_ref;
my @blast;

#declare hash
my %gi_hash;

use strict;
use warnings;

open (IN,"<",$ARGV[0]) || die "Error cannot read from gi.txt: $!\n";
@gi_ref = <IN>;
close IN;

foreach $gi_ref (@gi_ref) {
	chomp $gi_ref;
	$gi_hash{$gi_ref} = 1;
}

#test
while (my ($key, $value)  = each (%gi_hash)) {
	print "key: $key\t value: $value\n";
}

print "Enter name of directory to grab blastn files from (including final /):\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error cannot read from directory: $!\n";

open (MERGE,">>","filtered.blastn") || die "Error cannot write to merged.blastn: $!\n";

while ($file = readdir(DIR)) {
	chomp $file;
	#print "$file\n";#test
	if ($file =~ /\d+\.fasta\.blastn/) {
		$file =~ /(\d+)\.fasta\.blastn/;
		$gi = $1;
		#print "$gi\n";#test
		if ($gi_hash{$gi}) {
			print "found match\n";
			$path_to_file = $dir.$file;
			open (BLAST,"<",$path_to_file) || die "Error cannot read blastn file:
			$!\n";
			@blast = <BLAST>;
			close BLAST;
			
			while ($blast[$i]) {
				$line = $blast[$i];
				chomp $line;
				print MERGE "$line\n";
				$i++;
			}
			$i=0;	
		}
	}
}
close DIR;
close MERGE;
