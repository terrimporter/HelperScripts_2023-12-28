#!/usr/bin/perl
#April 1, 2014 edited to handle hundreds/thousands of FC OTUs and compare these to hundreds/thousands of BR OTUs
#Uses BLAST to compare FC OTUS against a database created from BR OTUs, then parses tabular output to figure out the pairings
#Still use SeqPrep to do the final pairings
#
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
my $FCdir; ### STDIN
my $BRdir; ### STDIN
my $originaldir;
my $cwd;
my $i=0;
my $file;
my $well;
my $output;
my $dbasename;
my $filesize;
my $FCpath;
my $BRpath;
my $blastfile;
my $dbasefile;
my $j=0;
my $line;
my $qseqid;
my $sseqid;
my $pident;
my $olength;
my $mismatch;
my $gapopen;
my $qstart;
my $qend;
my $sstart;
my $send;
my $evalue;
my $bitscore;
my $pident_cutoff = 90; ### edit BLAST match cutoffs here (lenient settings here) ###
my $olength_cutoff = 25;
my $mismatch_cutoff = 2;
my $gapopen_cutoff = 2;
my $well_id;
my $size;
my $percent;
my $id;
my $FCfile;
my $BRfile;
my $FC_header;
my $BR_header;
my $FC_seq;
my $BR_seq;
my $FC_fastq;
my $BR_fastq;
my $length;
my $outfile;
my $seqin;
my $seqout;
my $seqobj;
my $header;
my $seq;
my $FC_id;
my $BR_id;
my $FC_well;
my $BR_well;
my $FWD;
my $REV;
my $rc;
my $counter=0;
my $o = 80; ### Edit SeqPrep overlap paramaters here (stringent settings here) ###
my $n = 0.98;
my $m = 0.02;
my $newfile;

#declare array
my @FCfiles;
my @file;
my @BRfiles;
my @output;
my @line;
my @blastfile;
my @qseqid;
my @well_id;
my @sseqid;
my @blastfiles;

#declare hash
my %FC; #key = well; value = filename; append $FCdir path to file later
my %BR;
my %FC2;
my %BR2;
my %FC3;
my %BR3;
my %blastMatch; #key=qseqid; value=sseqid;

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
chomp $cwd;
$FCpath = $cwd;
#print "cwd: $cwd\n";#test
@FCfiles = qx(ls | grep .centroids.percent);

#change back to original dir
chdir $originaldir;
$cwd = cwd;
print "cwd: $cwd\n";

print "\nEnter path to directory containing *.centroids.percent files for one BR plate:\n";
$BRdir = <STDIN>;
chomp $BRdir;
print "\n";

#collect BR files and create a custom BLAST database for each well #NEW 04/01/14
chdir $BRdir;
$cwd = cwd;
chomp $cwd;
$BRpath = $cwd;
#print "cwd: $cwd\n";#test
@BRfiles = qx(ls | grep .centroids.percent);
#print "@BRfiles\n";#test

while ($BRfiles[$i]) {
	$file = $BRfiles[$i];
	chomp $file;
	$filesize = -s $file;

	if ($filesize > 0) {

		@file = split(/\./,$file);
		$well = $file[0];

		#create a custom BLAST 2.2.29 dbase for each BR file
		$dbasename = $well;
		$output = qx(makeblastdb -in $file -dbtype nucl -title $dbasename -out $dbasename);
	}

	$i++;
	$file=();
	@file=();
	$well=();
	$dbasename=();
	$output=();
	$filesize=();
}
$i=0;

#change back to original dir
chdir $originaldir;
$cwd = cwd;
print "cwd: $cwd\n";

chdir $FCdir;
$cwd = cwd;
print "\ncwd: $cwd\n";

#BLAST each FC file against corresponding BR dbase
while ($FCfiles[$i]) {
	$file = $FCfiles[$i];
	chomp $file;
	$filesize = -s $file;

	if ($filesize > 0 ) {

		@file = split(/\./,$file);
		$well = $file[0];
		$dbasename = $BRpath."/".$well;
		$dbasefile = $dbasename.".nsq";

		if (-e $dbasefile) {#check if corresponding BR dbase exists

			@blastfiles = BLAST_FC_vs_BRdbase(\$file,\$dbasename);
		}
	
	}

	$i++;
	$file=();
	@file=();
	$well=();
	$dbasename=();
	$dbasefile=();
	$filesize=();
	@blastfiles=();
}
$i=0;

parse_tabular_BLAST();

#test
while( ($qseqid,$sseqid) = each (%blastMatch) ) {
	print "FC: $qseqid\tBR: $sseqid\n";
}

#for each blast match, grab FC file, create part files as needed, grab BR file, create part files as needed, convert fasta to fastq, do seqprep pairing

while ( ($qseqid,$sseqid) = each(%blastMatch) ) {
	$counter++;

	#create FC part file
	@qseqid = split(/;/,$qseqid);
	$well_id = $qseqid[0];
	$size = $qseqid[1];
	$percent = $qseqid[2];

	@well_id = split(/_/,$well_id);
	$FC_well = $well_id[0];
#	$counter = $well_id[1];

	$FCfile = $FC_well.".fasta.trimmed.sorted.centroids.percent";
#	print "FCfile: $FCfile\n";#test - ok
	
	@output = qx(grep '$qseqid' $FCfile -A 1);	
	print "TEST:@output\n";#test - ok
	$FC_header = $output[0];
	chomp $FC_header;
	$FC_header =~ s/^>//;
	$FC_seq = $output[1];
	chomp $FC_seq;

	$outfile = $FCfile.".".$counter.".fasta"; #ok
	open (FCOUT, ">>", $outfile) || die "Error cannot open $outfile: $!\n";
	print FCOUT ">$FC_header\n$FC_seq\n";
	close FCOUT;

	#create BR part file
	@sseqid = split(/;/,$sseqid);
	$well_id = $sseqid[0];
	$size = $sseqid[1];
	$percent = $sseqid[2];

	@well_id = split(/_/,$well_id);
	$BR_well = $well_id[0];
#	$BR_id = $well_id[1];

	$BRfile = $BRpath."/".$BR_well.".fasta.trimmed.sorted.centroids.percent";
#	print "BRfile: $BRfile\n";#test

	@output = qx(grep '$sseqid' $BRfile -A 1);
#	print "@output\n";#test
	$header = $output[0];
	chomp $header;
	$seq = $output[1];
	chomp $seq;

	$outfile = $BRfile.".".$counter.".fasta";
	open (BROUT, ">>", $outfile) || die "Error cannto open $outfile: $!\n";
	print BROUT "$header\n$seq\n";
	close BROUT;

	#create BR fastq, needs to be reverse-complemented!
	$file = $BRpath."/".$BR_well.".fasta.trimmed.sorted.centroids.percent.".$counter.".fasta";
	$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$file);
	$seqout = Bio::SeqIO->new(-format=>'Fasta', -file=>'>>temp.fa');

	while ($seqobj = $seqin->next_seq() ) {
		$BR_header = $seqobj->display_id;
		$rc = $seqobj->revcom; #reverse complemented
		$seqout->write_seq($rc);
		$seqout=();
	}
	$seqobj=();
	$seqin=();

	open (TEMP, "<", 'temp.fa') || die "Error cannot open temp.fa: $!\n";

	while(<TEMP>){
		$line = $_;
		chomp $line;
		$length = length($line);

		if ($line !~ /^>/ && $length > 0) {
			push(@line,$line);
		}
	}
	$BR_seq = join('',@line);
	print "BR_seq: $BR_seq\n";
	close TEMP;

	#create FC fastq
	$FC_fastq = $originaldir."/FC.".$FCfile.".".$counter.".fastq";
	open (FCFASTQ, ">>", $FC_fastq) || die "Error cannot open $FC_fastq:$!\n";
	$length = length($FC_seq);
	print FCFASTQ "\@$FC_header|$BR_header\n";
	print FCFASTQ "$FC_seq\n";
	print FCFASTQ "+\n";
	print FCFASTQ 'I' x $length;
	print FCFASTQ "\n";
	close FCFASTQ;

	#create BR fastq
	$BR_fastq = $originaldir."/BR.".$BR_well.".fasta.trimmed.sorted.centroids.percent.".$counter.".fastq";
	open (BRFASTQ, ">>", $BR_fastq) || die "Error cannot open $BR_fastq:$!\n";
	$length = length($BR_seq);
	print BRFASTQ "\@$FC_header|$BR_header\n";
	print BRFASTQ "$BR_seq\n";
	print BRFASTQ "+\n";
	print BRFASTQ 'I' x $length;
	print BRFASTQ "\n";
	close BRFASTQ;

	unlink('temp.fa');
	@line=();

	chdir $originaldir;
	$cwd = cwd;
	print "\ncwd: $cwd\n";

	$FWD = "FC.".$FCfile.".".$counter.".fastq";
	$REV = "BR.".$BR_well.".fasta.trimmed.sorted.centroids.percent.".$counter.".fastq";
	run_SeqPrep(\$FWD,\$REV);

	chdir $FCdir;

}

####################

sub BLAST_FC_vs_BRdbase {

#declare var
my $file = ${$_[0]};
my $dbasename = ${$_[1]};
my $outfile;
my $output;
my $infile;
my $seqin;
my $seqout;
my $seqobj;
my $header;
my $indv_filename;
my $i=0;
my $seq;

#declare array
my @indv_files;
my @outfile;

$seqin = Bio::SeqIO->new(-format=>'Fasta', -file=>$file);

#split the original .percent file into an individual file for each OTU
while ($seqobj = $seqin->next_seq()) {
	$header = $seqobj->display_id;
	$newfile = $header;
	$newfile =~ s/;/_/;
	$newfile =~ s/;/_/;
	$newfile =~ s/;//;
	$newfile =~ s/=/_/g;
	$indv_filename = $newfile.".indv";

	open (INDV, ">>", $indv_filename) || die "Error cannot open $indv_filename:$!\n";

#	if (-e $indv_filename) {
#		$seqout = Bio::SeqIO->new(-format=>'Fasta', -file=> $indv_filename);
	push(@indv_files,$indv_filename);
	$seq=$seqobj->seq();
	print INDV ">$header\n$seq\n";
	close INDV;
#	}
#	else {
#		print "Error, could not open $indv_filename\n";
#	}

	$header=();
	$newfile=();
	$seq=();
#	$seqout=();
	$indv_filename=();
}

#BLAST each individual OTU against the database
while($indv_files[$i]) {
	$infile = $indv_files[$i];
	$outfile = $infile.".out";

	$output = system("blastn -query $infile -task megablast -db $dbasename -out $outfile -outfmt 6");

	push(@outfile, $outfile);
		
	$i++;
	$infile=();
	$outfile=();
	$output=();
}
$i=0;

return @outfile;

}

####################

sub parse_tabular_BLAST {

#parse blast *.out files, for good hits that matches overlap criteria
@output = qx(ls | grep .out);

while ($output[$i]) {
	$blastfile = $output[$i];
	chomp $blastfile;
	$filesize = -s $blastfile;

	if ( $filesize > 0 ) {
		open (BLAST, "<", $blastfile) || die "Error cannot open $blastfile\n";
		@blastfile = <BLAST>;
		close BLAST;

		while ($blastfile[$j]) {
			$line = $blastfile[$j];
			chomp $line;

			@line = split(/\t/,$line);
			$qseqid = $line[0];
			$sseqid = $line[1];
			$pident = $line[2];
			$olength = $line[3];
			$mismatch = $line[4];
			$gapopen = $line[5];
			$qstart = $line[6];
			$qend = $line[7];
			$sstart = $line[8];
			$send = $line[9];
			$evalue = $line[10];
			$bitscore = $line[11];

			if ($pident >= $pident_cutoff && $olength >= $olength_cutoff && $mismatch <= $mismatch_cutoff && $gapopen <= $gapopen_cutoff) {
				$blastMatch{$qseqid} = $sseqid;
			}
			
			$j++;
			$line=();
			@line=();
			$qseqid=();
			$sseqid=();
			$pident=();
			$olength=();
			$mismatch=();
			$gapopen=();
			$qstart=();
			$qend=();
			$sstart=();
			$send=();
			$evalue=();
			$bitscore=();

		}
		$j=0;
	}
	$i++;
	$blastfile=();
	$filesize=();
	@blastfile=();

}
$i=0;

}

####################

sub run_SeqPrep {

my $FWD = ${$_[0]};
my $REV = ${$_[1]};

my $FWDout = $FWD.".out";
my @FWDout = split(/\./,$FWDout); #parse to get frag, well, & id
my $frag = $FWDout[0];
my $well = $FWDout[1];
my $counter = $FWDout[7];
my $FWDprefix = $well.".".$frag.".".$counter;

my $REVout = $REV.".out";
my @REVout = split(/\./,$REVout);
$frag = $REVout[0];
$well = $REVout[1];
$counter = $REVout[7];
my $REVprefix = $frag.".".$counter;

my $PAIREDout = $FWDprefix.".".$REVprefix.".paired.fastq.gz";
my $PAIREDout2 = $FWDprefix.".".$REVprefix.".paired.aln.gz";
my $output;

$output = qx(SeqPrep -f $FWD -r $REV -1 $FWDout -2 $REVout -q 20 -s $PAIREDout -E $PAIREDout2 -o $o -n $n -m $m);

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
