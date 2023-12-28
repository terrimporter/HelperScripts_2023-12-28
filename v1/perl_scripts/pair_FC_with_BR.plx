#!/usr/bin/perl
#March 28, 2014 edit to work with 3 dominant otus now, 1:3*, 3:1, 3:3
##March 14, 2014 by Terri Porter
#Script to pair 99% sequence similarity delimited OTUs from FC and BR fragments
#For now, only work with wells that have a single domininant OTU (i.e. >20% reads per OTU)
#Add modules to handle FR:BR dominant OTU combinations of 1:2, 2:1, 2:2, etc. later
#USAGE perl pair_FR_with_BR.plx

use strict;
use warnings;
use Cwd;
use Bio::SeqIO;

#declare var
my $FCdir; ### user input required
my $BRdir; ### user input required
my $originaldir;
my $cwd;
my $i=0;
my $file;
my $numDomOTUs;
my $well;
my $FCfile;
my $BRfile;
my $file2;
my $fwd;
my $rev;
my $count=0;
my $FCfile1;
my $FCfile2;
my $FCfile3;
my $BRfile1;
my $BRfile2;
my $BRfile3;
my $FCversion;
my $BRversion;
my $seqin;
my $seqobj;
my $header;
my $seq;
my $flag=0; # 11 indicates special SeqPrep settings for FC:BR 1:1

#declare array
my @FCfiles;
my @file;
my @BRfiles;

#declare hash
my %FC; #key = well; value = filename; append $FCdir path to file later
my %BR;
my %FC2;
my %BR2;
my %FC3;
my %BR3;

$originaldir = cwd;
print "originaldir: $originaldir\n";#test

### Process one FC & BR plate at a time for now ###
print "\nEnter path to directory containing *.centroids.percent files for one FC plate:\n";
$FCdir = <STDIN>;
chomp $FCdir;
print "\n";

#collect FC files
chdir $FCdir;
$cwd = cwd;
#print "cwd: $cwd\n";#test

@FCfiles = qx(ls | grep .centroids.percent);
#print "@FCfiles\n";#test

while ($FCfiles[$i]) {
	$file = $FCfiles[$i];
	chomp $file;

	@file = split(/\./,$file);
	$well = $file[0];
	$numDomOTUs = qx(grep ">" $file | wc -l);

	### only deal with wells that contain 1, 2, or 3 dominant 99% OTU right now
	if ($numDomOTUs == 1 ) { 
		$FC{$well} = $file;
	}
	elsif ($numDomOTUs == 2) {
		$FC2{$well} = $file;
	}
	elsif ($numDomOTUs == 3) {
		$FC3{$well} = $file;
		print "added to hash FC3\n";#test
	}

	$i++;
	$file=();
	$numDomOTUs=();
	@file=();
	$well=();
}
$i=0;

#change back to original dir
chdir $originaldir;
$cwd = cwd;
print "cwd: $cwd\n";

print "\nEnter path to directory containing *.centroids.percent files for one BR plate:\n";
$BRdir = <STDIN>;
chomp $BRdir;
print "\n";

#collect BR files
chdir $BRdir;
$cwd = cwd;
#print "cwd: $cwd\n";#test

@BRfiles = qx(ls | grep .centroids.percent);
#print "@BRfiles\n";#test

while ($BRfiles[$i]) {
	$file = $BRfiles[$i];
	chomp $file;

	@file = split(/\./,$file);
	$well = $file[0];
	$numDomOTUs = qx(grep ">" $file | wc -l);

	### only deal with wells that contain 1, 2 or 3 dominant 99% OTU right now
	if ($numDomOTUs == 1 ) { 
		$BR{$well} = $file;
	}
	elsif ($numDomOTUs == 2) {
		$BR2{$well} = $file;
	}
	elsif ($numDomOTUs == 3) {
		$BR3{$well} = $file;
		print "added to hash BR3\n";#test
	}

	$i++;
	$file=();
	$numDomOTUs=();
	@file=();
	$well=();
}
$i=0;

#change back to original dir
chdir $originaldir;
$cwd = cwd;
print "cwd: $cwd\n";

### create all necessary part files first!!! ###
#don't worry about %FC
while ( ($well, $file) = each (%FC2) ) {
	$FCfile = $FCdir."/".$file;
	$FCfile1 = $FCfile.".1";
	$FCfile2 = $FCfile.".2";

	open (OUT1, ">>", $FCfile1) || die "Error cannot open $FCfile1: $!\n";
	open (OUT2, ">>", $FCfile2) || die "Error cannot open $FCfile2: $!\n";
	
	$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$FCfile);

	while ($seqobj = $seqin->next_seq()) {
		$count++;
		$header = $seqobj->display_id;

		$seq = $seqobj->seq;

		if ($count == 1) {
			print OUT1 ">$header\n$seq\n";
			close OUT1;
		}
		elsif ($count == 2) {
			print OUT2 ">$header\n$seq\n";
			close OUT2;
		}

		$header=(); 
		$seq=();
	}
	$count=0;
	$seqin=();

}

while ( ($well, $file) = each (%FC3) ) {
	$FCfile = $FCdir."/".$file;
	$FCfile1 = $FCfile.".1";
	$FCfile2 = $FCfile.".2";
	$FCfile3 = $FCfile.".3";

	open (OUT1, ">>", $FCfile1) || die "Error cannot open $FCfile1: $!\n";
	open (OUT2, ">>", $FCfile2) || die "Error cannot open $FCfile2: $!\n";
	open (OUT3, ">>", $FCfile3) || die "Error cannot open $FCfile3: $!\n";

	$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$FCfile);

	while ($seqobj = $seqin->next_seq()) {
		$count++;
		$header = $seqobj->display_id;

		$seq = $seqobj->seq;

		if ($count == 1) {
			print OUT1 ">$header\n$seq\n";
			close OUT1;
		}
		elsif ($count == 2) {
			print OUT2 ">$header\n$seq\n";
			close OUT2;
		}
		elsif ($count == 3) {
			print OUT3 ">$header\n$seq\n";
			close OUT3;
		}
		$header=(); 
		$seq=();
	}
	$count=0;
	$seqin=();

	close OUT1;
	close OUT2;
	close OUT3;
}
#don't worry about %BR
while ( ($well, $file) = each (%BR2) ) {
	$BRfile = $BRdir."/".$file;
	$BRfile1 = $BRfile.".1";
	$BRfile2 = $BRfile.".2";

	open (OUT1, ">>", $BRfile1) || die "Error cannot open $BRfile1: $!\n";
	open (OUT2, ">>", $BRfile2) || die "Error cannot open $BRfile2: $!\n";

	$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$BRfile);

	while ($seqobj = $seqin->next_seq()) {
		$count++;
		$header = $seqobj->display_id;

		$seq = $seqobj->seq;

		if ($count == 1) {
			print OUT1 ">$header\n$seq\n";
			close OUT1;
		}
		elsif ($count == 2) {
			print OUT2 ">$header\n$seq\n";
			close OUT2;
		}

		$header=(); 
		$seq=();
	}
	$count=0;
	$seqin=();

}

while ( ($well, $file) = each (%BR3) ) {
	$BRfile = $BRdir."/".$file;
	$BRfile1 = $BRfile.".1";
	$BRfile2 = $BRfile.".2";
	$BRfile3 = $BRfile.".3";

	open (OUT1, ">>", $BRfile1) || die "Error cannot open $BRfile1: $!\n";
	open (OUT2, ">>", $BRfile2) || die "Error cannot open $BRfile2: $!\n";
	open (OUT3, ">>", $BRfile3) || die "Error cannot open $BRfile3: $!\n";

	$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$BRfile);

	while ($seqobj = $seqin->next_seq()) {
		$count++;
		$header = $seqobj->display_id;

		$seq = $seqobj->seq;

		if ($count == 1) {
			print OUT1 ">$header\n$seq\n";
			close OUT1;
		}
		elsif ($count == 2) {
			print OUT2 ">$header\n$seq\n";
			close OUT2;
		}
		elsif ($count == 3) {
			print OUT3 ">$header\n$seq\n";
			close OUT3;
		}
		$header=(); 
		$seq=();
	}
	$count=0;
	$seqin=();
}

#only compare wells where FC has only 1 dominant OTU
while ( ($well,$file) = each (%FC) ) {
	if (exists $BR{$well}) { #FC:BR 1:1
		print "\nTrying to process FC:BR 1:1\n"; #test
		$file2 = $BR{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent
		$BRfile = $BRdir."/".$file2;
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile,\$BRfile,\$FCversion,\$BRversion); ### need headers to be same
		$flag=11;
#		print "fwd:$fwd\nrev:$rev\n"; #test
		run_SeqPrep(\$fwd,\$rev);
	
	}
	elsif (exists $BR2{$well}) { #FC:BR 1:2
		print "\nTrying to process FC:BR 1:2\n"; #test
		$file2 = $BR2{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent
		$BRfile = $BRdir."/".$file2; #split this file into two separate files
#		print "BRfile: $BRfile\n"; #test - ok

		$BRfile1 = $BRfile.".1";
		$BRfile2 = $BRfile.".2";

		#handle 1:file1 comparison
		print "Trying to handle 1:2a comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile,\$BRfile1,\$FCversion,\$BRversion);
#		print "fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle 1:file2 comparison
		print "Trying to handle 1:2b comparison now\n";
		$FCversion = 1;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile,\$BRfile2,\$FCversion,\$BRversion);
#		print "fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);
	}
	elsif (exists $BR3{$well}) { #FC:BR 1:3
		print "\nTrying to process FC:BR 1:3\n"; #test
		$file2 = $BR3{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent
		$BRfile = $BRdir."/".$file2; #split this file into two separate files
#		print "BRfile: $BRfile\n"; #test - ok

		$BRfile1 = $BRfile.".1";
		$BRfile2 = $BRfile.".2";
		$BRfile3 = $BRfile.".3";

		#handle 1:file1 comparison
		print "Trying to handle 1:2a comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile,\$BRfile1,\$FCversion,\$BRversion);
#		print "fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle 1:file2 comparison
		print "Trying to handle 1:2b comparison now\n";
		$FCversion = 1;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile,\$BRfile2,\$FCversion,\$BRversion);
#		print "fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle 1:file3 comparison
		print "Trying to handle 1:2c comparison now\n";
		$FCversion = 1;
		$BRversion = 3;
		($fwd,$rev) = fasta_to_fastq(\$FCfile,\$BRfile3,\$FCversion,\$BRversion);
#		print "fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);
	}

}

#now compare when FC has 2 dominant OTUs
while ( ($well,$file) = each (%FC2) ) {

	if (exists $BR{$well}) { #FC:BR 2:1
		print "\nTrying to process FC:BR 2:1\n"; #test
		$file2 = $BR{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent #split this file into two separate files
		$BRfile = $BRdir."/".$file2;

		$FCfile1 = $FCfile.".1";
		$FCfile2 = $FCfile.".2";

		#handle file1:1 comparison
		print "Trying to handle 1a:1 comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:1 comparison
		print "Trying to handle 1b:1 comparison now\n";
		$FCversion = 2;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

	}
	elsif (exists $BR2{$well}) { #FC:BR 2:2 ### problem here!!! ###
		print "\nTrying to process FC:BR 2:2\n"; #test
		$file2 = $BR2{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent # already split into 2 files
		print "FCfile: $FCfile\n";#test
		$BRfile = $BRdir."/".$file2; # already split into 2 files
#		print "BRfile: $BRfile\n";#test

		$FCfile1 = $FCfile.".1";
		if (-e $FCfile1) {
			print "FCfile1: $FCfile1 exists\n";
		}
		else {
			print "FCfile1: $FCfile1 does not exist\n";
		}
		$FCfile2 = $FCfile.".2";

		$BRfile1 = $BRfile.".1";
		$BRfile2 = $BRfile.".2";
		
		#handle file1:file1 comparison
		print "Trying to handle 1a:1a comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file1:file2 comparison
		print "Trying to handle 1a:1b comparison now\n";
		$FCversion = 1;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file1 comparison
		print "Trying to handle 1b:1a comparison now\n";
		$FCversion = 2;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file2 comparison
		print "Trying to handle 1b:2b comparison now\n";
		$FCversion = 2;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);
	}
	elsif (exists $BR3{$well}) { #FC:BR 2:3
		print "\nTrying to process FC:BR 2:2\n"; #test
		$file2 = $BR3{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent # already split into 2 files
		$BRfile = $BRdir."/".$file2; # already split into 3 files

		$FCfile1 = $FCfile.".1";
		$FCfile2 = $FCfile.".2";
		$BRfile1 = $BRfile.".1";
		$BRfile2 = $BRfile.".2";
		$BRfile3 = $BRfile.".3";
		
		#handle file1:file1 comparison
		print "Trying to handle 1a:1a comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file1:file2 comparison
		print "Trying to handle 1a:1b comparison now\n";
		$FCversion = 1;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file1:file3 comparison
		print "Trying to handle 1a:1c comparison now\n";
		$FCversion = 1;
		$BRversion = 3;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile3,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file1 comparison
		print "Trying to handle 1b:1a comparison now\n";
		$FCversion = 2;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file2 comparison
		print "Trying to handle 1b:2b comparison now\n";
		$FCversion = 2;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file3 comparison
		print "Trying to handle 1b:2b comparison now\n";
		$FCversion = 2;
		$BRversion = 3;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile3,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);
	}

}

#now compare when FC has 3 dominant OTUs
while ( ($well,$file) = each (%FC3) ) {

	if (exists $BR{$well}) { #FC:BR 3:1
		print "\nTrying to process FC:BR 3:1\n"; #test
		$file2 = $BR{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent #split this file into three separate files
		$BRfile = $BRdir."/".$file2;

		$FCfile1 = $FCfile.".1";
		$FCfile2 = $FCfile.".2";
		$FCfile3 = $FCfile.".3";

		#handle file1:1 comparison
		print "Trying to handle 1a:1 comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:1 comparison
		print "Trying to handle 1b:1 comparison now\n";
		$FCversion = 2;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file3:1 comparison
		print "Trying to handle 1c:1 comparison now\n";
		$FCversion = 3;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);
	}
	elsif (exists $BR2{$well}) { #FC:BR 3:2
		print "\nTrying to process FC:BR 3:2\n"; #test
		$file2 = $BR2{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent # already split into 2 files
		$BRfile = $BRdir."/".$file2; # already split into 2 files

		$FCfile1 = $FCfile.".1";
		$FCfile2 = $FCfile.".2";
		$FCfile3 = $FCfile.".3";
		$BRfile1 = $BRfile.".1";
		$BRfile2 = $BRfile.".2";
		
		#handle file1:file1 comparison
#		print "Trying to handle 1a:1a comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file1:file2 comparison
#		print "Trying to handle 1a:1b comparison now\n";
		$FCversion = 1;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file1 comparison
#		print "Trying to handle 1b:1a comparison now\n";
		$FCversion = 2;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file2 comparison
#		print "Trying to handle 1b:2b comparison now\n";
		$FCversion = 2;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file3:file1 comparison
#		print "Trying to handle 1c:1a comparison now\n";
		$FCversion = 3;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile3,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file3:file2 comparison
#		print "Trying to handle 1c:2b comparison now\n";
		$FCversion = 3;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile3,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

	}
	elsif (exists $BR3{$well}) { #FC:BR 3:3
		print "\nTrying to process FC:BR 3:3\n"; #test
		$file2 = $BR3{$well};
		$FCfile = $FCdir."/".$file; # *.centroids.percent # already split into 2 files
		$BRfile = $BRdir."/".$file2; # already split into 2 files

		$FCfile1 = $FCfile.".1";
		$FCfile2 = $FCfile.".2";
		$FCfile3 = $FCfile.".3";
		$BRfile1 = $BRfile.".1";
		$BRfile2 = $BRfile.".2";
		$BRfile3 = $BRfile.".3";
		
		#handle file1:file1 comparison
		print "Trying to handle 1a:1a comparison\n";
		$FCversion = 1;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file1:file2 comparison
		print "Trying to handle 1a:1b comparison now\n";
		$FCversion = 1;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file1:file3 comparison
		print "Trying to handle 1a:1c comparison now\n";
		$FCversion = 1;
		$BRversion = 3;
		($fwd,$rev) = fasta_to_fastq(\$FCfile1,\$BRfile3,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file1 comparison
		print "Trying to handle 1b:1a comparison now\n";
		$FCversion = 2;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file2 comparison
		print "Trying to handle 1b:2b comparison now\n";
		$FCversion = 2;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file2:file3 comparison
		print "Trying to handle 1b:2b comparison now\n";
		$FCversion = 2;
		$BRversion = 3;
		($fwd,$rev) = fasta_to_fastq(\$FCfile2,\$BRfile3,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file3:file1 comparison
		print "Trying to handle 1c:1a comparison now\n";
		$FCversion = 3;
		$BRversion = 1;
		($fwd,$rev) = fasta_to_fastq(\$FCfile3,\$BRfile1,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file3:file2 comparison
		print "Trying to handle 1c:2b comparison now\n";
		$FCversion = 3;
		$BRversion = 2;
		($fwd,$rev) = fasta_to_fastq(\$FCfile3,\$BRfile2,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);

		#handle file3:file3 comparison
		print "Trying to handle 1c:2b comparison now\n";
		$FCversion = 3;
		$BRversion = 3;
		($fwd,$rev) = fasta_to_fastq(\$FCfile3,\$BRfile3,\$FCversion,\$BRversion);
#		print "$fwd:$fwd\nrev:$rev\n"; #test
		$flag=0;
		run_SeqPrep(\$fwd,\$rev);
	}

}

clean_up_dir();

####################

sub clean_up_dir {

#declare var
my $file;

#declare array
my @output;

chdir $FCdir;
$cwd = cwd;
print "\nCleaning up $cwd\n";

@output = qx(ls | grep percent.1);
foreach $file (@output) {
	chomp $file;
	do {
		unlink $file and print "Deleted file: $file\n";
	} or warn "Cannot unlink file $file\n";
}
@output=();

@output = qx(ls | grep percent.2);
foreach $file (@output) {
	chomp $file;
	do {
		unlink $file and print "Deleted file: $file\n";
	} or warn "Cannot unlink file $file\n";
}
@output=();

@output = qx(ls | grep percent.3);
foreach $file (@output) {
	chomp $file;
	do {
		unlink $file and print "Deleted file: $file\n";
	} or warn "Cannot unlink file $file\n";
}
@output=();

chdir $originaldir;
$cwd = cwd;
#print "\n...back to original directory.\n";

chdir $BRdir;
$cwd = cwd;
print "\nNow cleaning up $cwd\n";

@output = qx(ls | grep percent.1);
foreach $file (@output) {
	chomp $file;
	do {
		unlink $file and print "Deleted file: $file\n";
	} or warn "Cannot unlink file $file\n";
}
@output=();

@output = qx(ls | grep percent.2);
foreach $file (@output) {
	chomp $file;
	do {
		unlink $file and print "Deleted file: $file\n";
	} or warn "Cannot unlink file $file\n";
}
@output=();

@output = qx(ls | grep percent.3);
foreach $file (@output) {
	chomp $file;
	do {
		unlink $file and print "Deleted file: $file\n";
	} or warn "Cannot unlink file $file\n";
}
@output=();

}

####################

sub run_SeqPrep {

#declare var
my $FWD = ${$_[0]};
my $REV = ${$_[1]};

my $FWDout = $FWD.".out"; #parse to get frag, version, & well
my @FWDout = split(/\./,$FWDout);
my $frag = $FWDout[0];
my $well = $FWDout[1];
my $version = $FWDout[2];
my $FWDprefix = $well.".".$frag.".".$version;

my $REVout = $REV.".out"; #parse to get frag & version
my @REVout = split(/\./,$REVout);
$frag = $REVout[0];
$well = $REVout[1];
$version = $REVout[2];
my $REVprefix = $frag.".".$version;

my $PAIREDout = $FWDprefix.".".$REVprefix.".paired.fastq.gz";
my $PAIREDout2 = $FWDprefix.".".$REVprefix.".paired.aln.gz";
my $output;

my $o = 25;
my $n = 1.0;
my $m = 0;

if ($flag == 11) { #FC:BR 1:1
	$o = 50;
	$n = 0.98;
	$m = 0.02;
}

$output = qx(SeqPrep -f $FWD -r $REV -1 $FWDout -2 $REVout -q 20 -s $PAIREDout -E $PAIREDout2 -o $o -n $n -m $m );

### options used for SeqPrep ###
# -f FC.fastq
# -r BR.fastq
# -1 FC.fastq.out (if trimming was done)
# -2 BR.fastq.out (if trimming was done)
# -q 20 (min Phred score of the overlap)
# -s paired.fastq.gz (the FC+BR paired sequence)
# -E paired.aln.gz (the FC+BR alignment)
# -o (minimum overlap length required CHANGE THIS TO SOMETHING HIGHER?!?)
# -n (default 0.9, min fraction of matching bases for overlap)
# -m (default 0.02, max number of good quality mismatches allowed in overlap)

$FWDout=();
$REVout=();
$PAIREDout=();
$PAIREDout2=();
$output=();

}

####################

sub fasta_to_fastq { ### right now, only works with one dominant taxon

use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;

#declare var
my $FWDin = ${$_[0]}; # *.centroids.percent
my $REVin = ${$_[1]};
my $FWDversion = ${$_[2]};
my $REVversion = ${$_[3]}; 
my $FWDout = "FC.".$well.".".$FWDversion.".fastq";
#print "FWDout:$FWDout\n"; #test -ok
my $REVout = "BR.".$well.".".$REVversion.".fastq";
#print "REVout: $REVout\n"; #test -ok

my $j;
my $line;
my $header;
my $seq;
my $length;
my $nextline;
my $rc; #revese-compliment
my $seqin;
my $seqobj;
my $seqout;
my $flag=0;
my $well_id;
my $size;
my $percent;
my $well;
my $id;
my $FC_frag;
my $BR_frag;
my $newHeader;

#declare array
my @rc;
my @seq;
my @temp;
my @header;
my @well_id;

#declare hash
my %FCheader; #key = well, value = FC_$id_$size_$percent
my %BRheader; #key = well, value = BR_$id_$size_percent

#process FC or forward reads first, for now just grab header data for later mapping
$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$FWDin);

while ($seqobj = $seqin->next_seq()) {
	$header = $seqobj->display_id;
	@header = split(/;/,$header);
	$well_id = $header[0];
	$size = $header[1];
	$percent = $header[2];

	$well_id =~ s/^>//g;
	@well_id = split(/_/,$well_id);
	$well = $well_id[0];
	$id = $well_id[1];

	$size =~ s/size=//g;
	$percent =~ s/percent=//g;

	$FCheader{$well} = "_FC_".$id."_".$size."_".$percent;
	
	$header=(); 
	@header=();
	$well_id=();
	$size=();
	$percent=();
	@well_id=();
	$well=();
	$id=();
}

$seqin=();

#process BR or reverse reads second, just to reverse-complement!
$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$REVin);
$seqout = Bio::SeqIO->new(-format=>'Fasta', -file=>'>>temp.fa');

while ($seqobj = $seqin->next_seq()) {
	$header = $seqobj->display_id;

	$rc=$seqobj->revcom;
	$seqout->write_seq($rc);

	$header=();
	$rc=();
	$seqout=();
}
$seqin=();
$seqout=();

#parse $REVin file, just to retrieve header data for later mapping
#$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$REVin);

#while ($seqobj = $seqin->next_seq()) {
#	$header = $seqobj->display_id;
#	@header = split(/;/,$header);
#	$well_id = $header[0];
#	$size = $header[1];
#	$percent = $header[2];

#	$well_id =~ s/^>//g;
#	@well_id = split(/_/,$well_id);
#	$well = $well_id[0];
#	$id = $well_id[1];

#	$size =~ s/size=//g;
#	$percent =~ s/percent=//g;

#	$BRheader{$well} = "_BR_".$id."_".$size."_".$percent;
	
#	$header=(); 
#	@header=();
#	$well_id=();
#	$size=();
#	$percent=();
#	@well_id=();
#	$well=();
#	$id=();
#}

#$seqin=();

#parse the temp.fa file, just to retrieve header data for later mapping
open (TEMP, "<", "temp.fa") || die "Error cannot open temp.fa: $!\n";
@temp = <TEMP>;
close TEMP;

while ($temp[$i]) {
	$line = $temp[$i];
	chomp $line;

	if ( $flag==0 && $line =~ /^>/) {
		$header = $line;
		@header = split(/;/,$header);
		$well_id = $header[0];
		$size = $header[1];
		$percent = $header[2];

		$well_id =~ s/^>//g;
		@well_id = split(/_/,$well_id);
		$well = $well_id[0];
		$id = $well_id[1];

		$size =~ s/size=//g;
		$percent =~ s/percent=//g;

		$BRheader{$well} = "_BR_".$id."_".$size."_".$percent;

	}
	$i++;

	$line=();
	$header=();
	@header=();
	$well_id=();
	$size=();
	$percent=();
	@well_id=();
	$well=();
	$id=();
}
$i=0;

#now use mapping files to create new merged headers for FC file first
if (-e $FWDout) {
	unlink $FWDout; 
	#trash old fastq with old merged header before creating new one!
}

open (OUT1, ">>", $FWDout) || die "Error cannot open $FWDout: $!\n";

$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$FWDin);

while ($seqobj = $seqin->next_seq()) {
	$header = $seqobj->display_id;
	@header = split(/;/,$header);
	$well_id = $header[0];

	$well_id =~ s/^>//g;
	@well_id = split(/_/,$well_id);
	$well = $well_id[0];

	if (exists $FCheader{$well}) {
		$FC_frag = $FCheader{$well};
		if (exists $BRheader{$well}) {
			$BR_frag = $BRheader{$well};
			$newHeader = $well.$FC_frag.$BR_frag;
		}
		else {
			print "Error processing BR well $well to create new merged header\n";
		}
	}
	else {
		print "Error processing FC well $well to create new merged header\n";
	}

	$seq = $seqobj->seq;
	@seq = split(//,$seq);
	$length = scalar(@seq);

	print OUT1 "\@$newHeader\n"; #tested, SeqPrep does not require illumina header format, FC and BR headers must match
	print OUT1 "$seq\n";
	print OUT1 "+\n";
	print OUT1 'I' x $length;
	print OUT1 "\n";

	$header=(); 
	@header=();
	$well_id=();
	@well_id=();
	$well=();
	$FC_frag=();
	$BR_frag=();
	$newHeader=();
	$seq=();
	@seq=();
	$length=();
}
close OUT1;
$seqin=();

#now use mapping files to create new merged headers for BR file second
#open (OUT2, ">>", $REVout) || die "Error cannot open $REVout: $!\n";

#$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$REVin);

#while ($seqobj = $seqin->next_seq()) {
#	$header = $seqobj->display_id;
#	@header = split(/;/,$header);
#	$well_id = $header[0];

#	$well_id =~ s/^>//g;
#	@well_id = split(/_/,$well_id);
#	$well = $well_id[0];

#	if (exists $FCheader{$well}) {
#		$FC_frag = $FCheader{$well};
#		if (exists $BRheader{$well}) {
#			$BR_frag = $BRheader{$well};
#			$newHeader = $well.$FC_frag.$BR_frag;
#		}
#		else {
#			print "Error processing BR well $well to create new merged header\n";
#		}
#	}
#	else {
#		print "Error processing FC well $well to create new merged header\n";
#	}

#	$seq = $seqobj->seq;
#	@seq = split(//,$seq);
#	$length = scalar(@seq);

#	print OUT2 "\@$newHeader\n"; #tested, SeqPrep does not require illumina header format, FC and BR headers must match
#	print OUT2 "$seq\n";
#	print OUT2 "+\n";
#	print OUT2 'I' x $length;
#	print OUT2 "\n";

#	$header=(); 
#	@header=();
#	$well_id=();
#	@well_id=();
#	$well=();
#	$FC_frag=();
#	$BR_frag=();
#	$newHeader=();
#	$seq=();
#	@seq=();
#	$length=();
#}
#close OUT2;
#$seqin=();

#now use mapping files to create new merged headers for BR fastq file (rc)
open (TEMP, "<", "temp.fa") || die "Error cannot open temp.fa: $!\n";
@temp = <TEMP>;
close TEMP;

if (-e $REVout) {
	unlink $REVout;
	#trash old fastq with old merged header before creating new one
}

open (OUT2, ">>", $REVout) || die "Error cannot open $REVout: $!\n";

while ($temp[$i]) {
	$line = $temp[$i];
	chomp $line;

	if ( $flag==0 && $line =~ /^>/) {
		$header = $line;
		@header = split(/;/,$header);
		$well_id = $header[0];

		$well_id =~ s/^>//g;
		@well_id = split(/_/,$well_id);
		$well = $well_id[0];

		if (exists $FCheader{$well}) {
			$FC_frag = $FCheader{$well};
			if (exists $BRheader{$well}) {
				$BR_frag = $BRheader{$well};
				$newHeader = $well.$FC_frag.$BR_frag;
			}
			else {
				print "Error processing BR well $well to create new merged header\n";
			}
		}
		else {
			print "Error processing FC well $well to create new merged header\n";
		}

		$flag=1;	
		print OUT2 "\@$newHeader\n";
	}

	elsif ( $flag==1 && $line !~ /^>/ ) {
		$seq = $line;
		$flag=2;
	}
	elsif ($flag==2 && $line !~ /^>/) {
		$seq = $seq.$line;
		$flag=2;
	}
	elsif ($flag==2 && $line =~ /^>/) {	
		print OUT2 "$seq\n";
		print OUT2 "+\n";
		@seq = split(//,$seq);
		$length = scalar(@seq);
		print OUT2 'I' x $length;
		print OUT2 "\n";
		$i--;
		$flag=0;
		@seq=();
		$seq=();
		$length=();
	}
	$i++;

	$line=();
	$header=();
	@header=();
	$well_id=();
	@well_id=();
	$well=();
	$FC_frag=();
	$BR_frag=();
	$newHeader=();
}
$i=0;

#don't forget to process the last sequence! (rc)
print OUT2 "$seq\n";
print OUT2 "+\n";
@seq = split(//,$seq);
$length=scalar(@seq);
print OUT2 'I' x $length;
print OUT2 "\n";
$flag=0;
@seq=();
$seq=();
$length=();
close OUT2;

unlink "temp.fa";

return ($FWDout, $REVout);

}

####################


