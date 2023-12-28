#!/usr/bin/perl
#Oct.20, 2011 by Terri Porter
#Script to grab gi from >gi\nseq\n
#usage perl get_gi_from_fasta_header.plx file.fasta

#declare var
my $i=0;
my $line;
my $gi;

#declare array
my @in;
my @gi;

open (IN,"<",$ARGV[0]) || die "Error cannot read file.fasta: $!\n";
@in = <IN>;
close IN;

open (OUT,">>","gi.txt") || die "Error cannot write to gi.txt: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>\d+/) {
		$line =~ /^>(\d+)/;
		$gi = $1;
		#push (@gi, $gi);
		print OUT "$gi\n";
	}
	$i++;
}
$i=0;
close OUT;
