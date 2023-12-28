#!/usr/bin/perl
#Nov. 1, 2012 by Terri Porter
#Script to find partial or missing *.out files so that they can be put in their own directory for another round of qsub_blast
#usage find_missing_files.plx

use warnings;
use strict;

#declare var
my $i=0;
my $filename;
my $prefix;
my $fileID;
my $value;
my $dirname;### edit this path ###
my $dirname2;### edit this path ###
my $outfile;
my $blastfile;
my $output;
my $batchsize=5000;### edit this number ###

#declare arry
my @outFiles;
my @outBlastnFiles;
my @filename;

#declare hash
my %outFiles;
my %outBlastnFiles;

#grab original query (*.out) files and hash them
@outFiles = qx(ls | grep '.out\$');
while ($outFiles[$i]) {
	$filename = $outFiles[$i];
	@filename = split(/\./, $filename);
	$prefix = $filename[0];
	$prefix =~ s/part//;
	$fileID = $prefix;
	$outFiles{$fileID} = 1;
	$i++;
	$filename=();
	@filename=();
	$prefix=();
	$fileID=();
}
$i=0;

#grab resulting blastn files (*.blastn) and hash them
@outBlastnFiles = qx(ls | grep '.out.blastn\$');
while ($outBlastnFiles[$i]) {
	$filename = $outBlastnFiles[$i];
	@filename = split(/\./, $filename);
	$prefix = $filename[0];
	$prefix =~ s/part//;
	$fileID = $prefix;
	$outBlastnFiles{$fileID} = 1;
	$i++;
	$filename=();
	@filename=();
	$prefix=();
	$fileID=();
}
$i=0;

#create new directory to hold *.out files that need to be re-blasted
$dirname = "/net/info254/1/home/terri/yersinia/R1_trimmed_MIN_MAX/all_blastn/redo";
mkdir $dirname, 0755;

#create new directory to hold *.out.blastn files that need to be deleted
$dirname2 = "/net/info254/1/home/terri/yersinia/R1_trimmed_MIN_MAX/all_blastn/delete";
mkdir $dirname2, 0755;

#for each outfile, check for presence and completeness of blastfile, mv outfiles that need to be re-blasted or blastfiles that need to be deleted
while (($fileID,$value) = each (%outFiles)) {
	$blastfile = "part".$fileID.".out.blastn";
	$outfile = "part".$fileID.".out";

	if (exists $outBlastnFiles{$fileID}) {
		$output = qx(grep "Query=" $blastfile | wc -l); #should be 25000 or 5000 if complete

		if ($output < $batchsize) {
			qx(mv $outfile $dirname);
			qx(mv $blastfile $dirname2);
		}
		elsif ($output > $batchsize) { # add this just in case 
			qx(mv $outfile $dirname);
			qx(mv $blastfile $dirname2);
		}
	}
	else {
		qx(mv $outfile $dirname);
	}
}
