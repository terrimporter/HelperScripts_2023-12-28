#!/usr/bin/perl
#Nov. 26, 2011 by Terri Porter
#Script to edit .qual (seqtrim.pl) for each seq in ITS.filtered (filter_fasta_by_length.plx) according to trimming info in outdata.csv (FungalITSextractor.pl)
#usage perl trim_qual_by_outdata.plx ITS.filtered outdata.csv seqtrim.qual.newlineremoved

#declare var
my $i=0;
my $line;
my $header;
my $id;
my $ITS1_trim;
my $start="nil";
my $end="nil";
my $id_qual;
my $j;
my $seq;
my $tot_length;
my $start_adj;
my $end_adj;
my $k=0;
my $base;

#declare array
my @fasta;
my @outdata;
my @qual;
my @ids;
my @line;
my @seq;

#declare hash
my %id_start;
my %id_end;

open (IN,"<", $ARGV[0]) || die "Error cannot open ITS.filtered: $!\n";
@fasta = <IN>;
close IN;

open (IN,"<",$ARGV[1]) || die "Error cannot open outdata.csv: $!\n";
@outdata = <IN>;
close IN;

open (IN,"<",$ARGV[2]) || die "Error cannot open seqtrim.qual: $!\n";
@qual = <IN>;
close IN;

#grab list of ids from ITS.filtered
while ($fasta[$i]) {
	$line = $fasta[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$id = substr $header, 1, 14;
		push (@ids,$id);
	}
	$i++;
}
$i=0;

#parse trim data from the whole outdata file
while ($outdata[$i]) {
	$line = $outdata[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$id = $line[0];
	$id =~ /\(.+\)\s+(\w+)/;
	$id = $1;
	$id = substr $id, 0, 14;	
	$ITS1_trim = $line[5];
	if ($ITS1_trim !~ /-----/) {
		$ITS1_trim =~ /ITS1:\s+(\d+)-(end|\d+)/;
		$start = $1;
		$end = $2;
		#print "start: $start\t end: $end\n";#test
	}
	else {
		$start = "nil";
		$end = "nil";
	}
	$id_start{$id} = $start;
	$id_end{$id} = $end;
	#test
	#print "id: $id\t trim_start: $start\t trim_end: $end\n";
	$i++;
}
$i=0;

open (OUT,">>","seqtrim.qual.newlineremoved.trimmed") || die "Error cannot write to seqtrim.qual.newlineremoved.trimmed: $!\n";

#for each fasta id, grab trim info from hashes, edit seqtrim.qual
while ($qual[$i]) {
	$line = $qual[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$header = $line;
		$header =~ /^>(\w+)\s+/;
		$id_qual = $1;
		#print "id_qual: $id_qual\n";#test
		$j = $i+1;
		$seq = $qual[$j];
		chomp $seq;
		#print "seq: $seq\n";#test
		@seq = split(/\s+/,$seq);
		$tot_length = scalar(@seq);
		#print "tot_length: $tot_length\n";#test
		#print "seq: @seq\n";#test

		$start = $id_start{$id_qual};
		#print "start: $start\n";#test
		if ($start eq "nil") {
			$i++;
			next;
		}
		else {
			$start_adj = $start-1;
			#print "start_adj: $start_adj\n";#test
		}

		$end = $id_end{$id_qual};
		if ($end eq "end") {
			$end = $tot_length;
			$end_adj = $end-1;
		}
		elsif ($end eq "nil") {
			$i++;
			next;
		}
		else {	
			$end_adj = $end-1;
		}

		while ($seq[$k]) {
			$base = $seq[$k];
			#print "base: $base\t k: $k\t start_adj: $start_adj\t end_adj: $end_adj\n";
			#print "$base ";
			if ($k >= $start_adj) {
				#print "$k >= $start_adj\n";
				if ($k <= $end_adj) {
					#print "$k<= $end_adj\n";
					push(@trimmed_seq,$base);
					#print "@trimmed_seq\n";
				}
				else {
					$k++;
					next;
				}
			}
			else {
				$k++;
				next;
			}
			$k++;
		}
		#print "@trimmed_seq\n";#test
		#print "\n";
		$k=0;
		$trimmed_seq = join(" ",@trimmed_seq);
		#print "trimmed_seq: $trimmed_seq\n";#test
		print OUT ">$id_qual\n$trimmed_seq\n";
	}
	$i++;
	@trimmed_seq=();
	$id_qual=();
	$seq=();
	$start="nil";
	$end="nil";
	$base=();
	$trimmed_seq=();
}
$i=0;
close OUT;
