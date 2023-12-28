#!/usr/bin/perl
#Nov.18,2011 by Terri Porter
#Script to organize a directory of prefix.sorted.fasta.trimmed and prefix.sorted.qual.trimmed into their own subdirectories called prefix
#usage perl directory_organize.plx

use strict;
use warnings;
use File::Copy;

#declare var
my $dir;
my $file;
my $prefix;
my $i=0;
my $path_to_newdir;
my $old_pathfile;
my $new_pathfile;
my $qualfile;

#declare array
my @filenames;
my @file;

#declare hash
my %prefix_fasta;
my %prefix_qual;

print "Enter path to directory containing files to be organized (including final /):\n";
$dir = <STDIN>;
chomp $dir;

#organize filenames in hashes indexed by prefix (sample name)
opendir (DIR, $dir) || die "Cannot open directory: $!\n";

while ($file = readdir(DIR)) {
	next if ($file =~ /^\./); #skip . and .. files
	@file = split(/\./,$file);
	$prefix = $file[0];
	if ($file =~ /\.fasta\.trimmed$/) {
		$prefix_fasta{$file} = $prefix;
	}
	elsif ($file =~ /\.qual\.trimmed$/) {
		$prefix_qual{$prefix} = $file;
	}
}

closedir DIR;

while (my ($file, $prefix) = each (%prefix_fasta)) {

	$path_to_newdir = $dir.$prefix."/";
	mkdir($path_to_newdir, 0777) || "Error canot create new dir: $!\n";
	
	$old_pathfile = $dir.$file; #move .fasta to new dir
	$new_pathfile = $path_to_newdir.$file;
	move($old_pathfile,$new_pathfile);
	
	$qualfile = $prefix_qual{$prefix}; #move .qual to new dir
	$old_pathfile = $dir.$qualfile;
	$new_pathfile = $path_to_newdir.$qualfile;
	move($old_pathfile,$new_pathfile);

}
