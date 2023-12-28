#!/usr/bin/perl
#TErri Porter, Feb.4.2011
#retrieve raw data records fron genbank, save, parse
#usage eutilities_gb2.plx gi.list

use strict;
use warnings;

use Bio::DB::EUtilities;
use Bio::SeqIO;
use Statistics::Lite qw(:all);

#var
my $factory;
my $filename;
my $seqin;
my $seq;
my $id;
my $feat_object;
my $sequence_fragment;
my $value;
my $SSU_seq="nil";
my $line;
my $i=0;
my $length;
my $min_length;
my $max_length;
my $mean_length;
my $id_from_record;

#array
my @ids;
my @id;
my @file2;
my @line;
my @SSU_seq;
my @length;

open (IN,"<",$ARGV[0]) || die ("Error cannot read gi.list: $!\n");
@ids=<IN>;
close IN;
foreach $id (@ids) {
	chomp $id;
}

open (OUT,">>","features.txt") || die ("Error cannot write to features.txt: $!\n");
#print OUT "ID\t18S sequence\n";

while ($ids[$i]) {
	$id = $ids[$i];
	push(@id,$id);
	$factory = Bio::DB::EUtilities -> new (-eutil => 'efetch',
						-db => 'nucleotide',
						-rettype => 'gb',
						-id => \@id);

	$filename = $id.".gb";

	$factory -> get_Response(-file => $filename);

	$seqin =Bio::SeqIO -> new(-file => $filename,
			-format => 'genbank');

	while ($seq = $seqin -> next_seq) {
#		$id_from_record = $seq->id;
		for $feat_object ($seq->get_SeqFeatures) {
			if ($feat_object->primary_tag eq "rRNA") {
				$sequence_fragment =  $feat_object->spliced_seq->seq;
				if ($feat_object->has_tag('product')) {
					for $value ($feat_object->get_tag_values('product')) {
						if ($value =~ /18S/) {
							$SSU_seq = $sequence_fragment;
						}
					}
				}
			}
		}
	}
	print OUT "$id\t$SSU_seq\n";
	$SSU_seq="nil";
	@id=();
	$i++;
	unlink($filename);
}
close OUT;

#parse features.txt and produce SSU fasta file, no "nil"

open (IN2,"<","features.txt") || die ("Error cannot read from features.txt:$!\n");
@file2 = <IN2>;
close IN2;

open (OUT2,">>","SSU_strict.fasta") || die ("Error cannot write to SSU_strict.fasta: $!\n");
$i=0;
while ($file2[$i]) {
	$line = $file2[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$id = $line[0];
	$SSU_seq = $line[1];
	if ($SSU_seq =~ /nil/) {
		$i++;
		next;
	}
	else {
		@SSU_seq = split(//,$SSU_seq);
		$length = scalar(@SSU_seq);
		push(@length,$length);
		print OUT2 ">$id|$length\n$SSU_seq\n";
	}
	$i++;
}
close OUT2;

$min_length = min(@length);
$max_length = max(@length);
$mean_length = mean(@length);

open (OUT3,">>","SSU_strict.stats") || die ("Error cannot write to SSU_strict.stats: $!\n");
print OUT3 "SSU stats for all sequences in Genbank with annotated 18S regions\n\n";
print OUT3 "Feature\tMin\tMax\tMean\n";
print OUT3 "SSU:\t$min_length\t$max_length\t$mean_length\n";
close OUT3;
