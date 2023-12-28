#!/usr/bin/perl
#Dec. 22, 2011 by Terri Porter
#Script to combine a bunch of phylip files from nex_to_phylip.plx into a single concatenated phylip file for prottest and raxml
#usage perl combine_phylip.plx

use strict;
use warnings;

#declare var
my $path;
my $i=0;
my $file;
my $phylip;
my $path_to_file;
my $j=0;
my $taxa;
my $char;
my $id;
my $seq;
my $sequence;
my $line;
my $new_sequence;
my $new_sum;
my $previous_sum;
my $path_to_outfile;

#declare array
my @files;
my @phylip;
my @in;
my @line;
my @line2;
my @char;

#declare hash
my %alignment;
my %concatenated;

print "Enter path to directory containing .phy files including final / : \n";
$path = <STDIN>;
chomp $path;

opendir (DIR,$path);
@files = readdir(DIR);
closedir (DIR);

#test
#print "files: @files\n";

while ($files[$i]) {
	$file = $files[$i];

	if ($file =~ /phy$/) {
		push(@phylip,$file);
	}
	$i++;
}
$i=0;

#test
#print "phylip files: @phylip\n";

while ($phylip[$i]) {
	$phylip = $phylip[$i];
	$path_to_file = $path.$phylip;
	
	open (IN,"<",$path_to_file) || die "Error cannot read phylip file: $!\n";
	@in = <IN>;
	close IN;

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($j==0) {
			@line = split(/\t/,$line);
			$taxa = $line[0];
			$char = $line[1];
			push(@char,$char);
			@line=();
			#print "taxa: $taxa characters: $char\n"; #test
		}
		else {
			@line2=split(/\s{8}/,$line);
			$id = $line2[0];
			$seq = $line2[1];
			$alignment{$id} = $seq;
			@line2=();
			#print "id: $id\t seq: $seq\n";#test
		}
		$j++;
	}
	$j=0;

	if ($i==0) {
		while (my ($id,$seq) = each (%alignment)) {
			$concatenated{$id} = $seq;
			#print "$id\t$seq\n";
		}
	}
	elsif ($i>=0) {
		while (my ($id,$seq) = each (%alignment)) {
			$sequence = $concatenated{$id};
			$new_sequence = $sequence.$seq;
			$concatenated{$id} = $new_sequence;
		}
	}

	$i++;
	%alignment=();
}
$i=0;

while($char[$i]) {
	$char = $char[$i];
	if ($i==0) {
		$previous_sum = $char;
	}
	else {
		$new_sum = $previous_sum+$char;
		$previous_sum = $new_sum;
	}
	$i++;
}
$i=0;

$path_to_outfile = $path."concatenated.phy";
open (OUT,">>",$path_to_outfile) || die "Error cannot write to outfile: $!\n";

#print "total char: $new_sum\n";
print OUT "$taxa\t$new_sum\n";
#test
while (my($id,$seq) = each (%concatenated)) {
	print OUT "$id        $seq\n";
}
close OUT;
