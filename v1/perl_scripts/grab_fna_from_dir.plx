#!/usr/bin/perl
#July 30, 2012 by Terri Porter
#Script to remove *.fna files from a list of directories and print them into a single dir
#usage perl grab_fna_from_dir.plx

use strict;
use warnings;

#declare var
my $topDirName;
my $i=0;
my $subDir;
my $pathToSubDir;
my $j=0;
my $subDirFile;
my $pathToSubDirFile;
my $k=0;
my $line;

#declare array
my @directoryList;
my @subDirList;
my @file;

print "Enter path to directory that contains the list of directories containing *.fna files (including final /:\n";

$topDirName = <STDIN>;
chomp $topDirName;

opendir (DIR, $topDirName) || die "Error in opening directory $topDirName\n";
@directoryList = readdir(DIR);
closedir(DIR);

open (OUT,">>", "merged.txt") || die "Error opening outfile: $!\n";

while ($directoryList[$i]) {
	$subDir = $directoryList[$i];
	chomp $subDir;

	if ($subDir !~ /\.$/) {
		$pathToSubDir = $topDirName.$subDir."/";
		opendir (SUBDIR, $pathToSubDir) || die "Error in opening subdirectory $pathToSubDir\n";
		@subDirList = readdir(SUBDIR);
		closedir(SUBDIR);
	}

	while ($subDirList[$j]) {
		$subDirFile = $subDirList[$j];
		chomp $subDirFile;
	
		if ($subDirFile !~ /\.$/) {
			#print "$subDirFile\n";
			if ($subDirFile =~ /fna$/) {
				$pathToSubDirFile = $pathToSubDir.$subDirFile;
				print "$pathToSubDirFile\n";
				open (IN, "<", $pathToSubDirFile) || die "Error in opening fna file: $!\n";
				@file = <IN>;
		
				while ($file[$k]) {
					$line = $file[$k];
					chomp $line;
					print OUT "$line\n";
					$k++;
				}
				$k=0;
				print OUT "\n";
				$line=();
			}
		}
		$j++;
		$subDirFile=();
		$pathToSubDirFile=();
		@file=();
	}
	$j=0;
	$i++;
	$subDir=();
	$pathToSubDir=();
	@subDirList=();
}
$i=0;
close OUT;
