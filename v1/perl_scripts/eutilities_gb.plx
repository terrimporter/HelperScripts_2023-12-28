#!/usr/bin/perl
#TErri Porter, Feb.4.2011
#retrieve raw data records fron genbank, save, parse

use strict;
use warnings;

use Bio::DB::EUtilities;
use Bio::SeqIO;
use Statistics::Lite qw(:all);

#var
my $factory;
my $file;
my $seqin;
my $seq;
my $id;
my $feat_object;
my $sequence_fragment;
my $value;
my $SSU_seq="nil";
my $ITS_seq="nil";
my $LSU_seq="nil";
my $spacer1="nil";
my $spacer2="nil";
my $line;
my $i=0;
my $length1;
my $length2;
my $length_ITS;
my $min_length1;
my $min_length2;
my $min_length_ITS;
my $max_length1;
my $max_length2;
my $max_length_ITS;
my $mean_length1;
my $mean_length2;
my $mean_length_ITS;
my $length1b;
my $length1b_;
my $length2b;
my $length2b_;
my $length_ITSb;
my $min_length1b;
my $min_length_ITSb;
my $min_length2b;
my $max_length1b;
my $max_length_ITSb;
my $max_length2b;
my $mean_length1b;
my $mean_length_ITSb;
my $mean_length2b;

#array
my @ids;
my @file2;
my @line;
my @spacer1;
my @spacer2;
my @ITS_seq;
my @length1;
my @length2;
my @length_ITS;
my @spacer1b;
my @spacer2b;
my @ITS_seqb;
my @length1b;
my @length2b;
my @length_ITSb;

open (IN,"<","gi.list") || die ("Error cannot read gi.list: $!\n");
@ids=<IN>;
close IN;
foreach $id (@ids) {
	chomp $id;
}


$factory = Bio::DB::EUtilities -> new (-eutil => 'efetch',
					-db => 'nucleotide',
					-rettype => 'gb',
					-id => \@ids);

$file = 'myseqs.gb';

$factory -> get_Response(-file => $file);

$seqin =Bio::SeqIO -> new(-file => $file,
			-format => 'genbank');

open (OUT,">>","features.txt") || die ("Error cannot write to features.txt: $!\n");
print OUT "ID\t18S\tspacer1\t5.8S\tspacer2\t28S\n";

while ($seq = $seqin -> next_seq) {
	$id = $seq->id;
	for $feat_object ($seq->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "rRNA") {
			$sequence_fragment =  $feat_object->spliced_seq->seq;
			if ($feat_object->has_tag('product')) {
				for $value ($feat_object->get_tag_values('product')) {
					if ($value =~ /18S/) {
						$SSU_seq = $sequence_fragment;
					}
					elsif ($value =~ /5.8S/) {
						$ITS_seq = $sequence_fragment;
					}
					elsif ($value =~ /(28S|26S|25S|large)/) {
						$LSU_seq = $sequence_fragment;
					}
				}
			}
		}
		elsif ($feat_object->primary_tag eq "misc_RNA") {
			$sequence_fragment =  $feat_object->spliced_seq->seq;
			if ($feat_object->has_tag('product')) {
				for $value ($feat_object->get_tag_values('product')) {
					if ($value =~ /spacer 1/) {
						$spacer1 = $sequence_fragment;
					}
					elsif ($value =~ /spacer 2/) {
						$spacer2 = $sequence_fragment;
					}
				}
			}
		}
	}
	print OUT "$id\t$SSU_seq\t$spacer1\t$ITS_seq\t$spacer2\t$LSU_seq\n";
	$SSU_seq="nil";
	$ITS_seq="nil";
	$LSU_seq="nil";
	$spacer1="nil";
	$spacer2="nil";
}
close OUT;

#parse features.txt and produce ITS only file, no "nil", count ITS1 and ITS2 >=100bp

open (IN2,"<","features.txt") || die ("Error cannot read from features.txt:$!\n");
@file2 = <IN2>;
close IN2;

open (OUT2,">>","ITS.fasta") || die ("Error cannot write to ITS.fasta: $!\n");

while ($file2[$i]) {
	$line = $file2[$i];
	chomp $line;
	@line = split(/\t/,$line);
	$id = $line[0];
	$SSU_seq = $line[1];
	$spacer1 = $line[2];
	$ITS_seq = $line[3];
	$spacer2 = $line[4];
	$LSU_seq = $line[5];
	if ($spacer1 =~ /nil/) {
		$i++;
		next;
	}
	elsif ($spacer2 =~ /nil/) {
		$i++;
		next;
	}
	elsif ($ITS_seq =~/nil/) {
		$i++;
		next;
	}
	else {
		@spacer1 = split(//,$spacer1);
		$length1 = scalar(@spacer1);
		push(@length1,$length1);
		@spacer2 = split(//,$spacer2);
		$length2 = scalar(@spacer2);
		push(@length2,$length2);
		@ITS_seq = split(//,$ITS_seq);
		$length_ITS = scalar(@ITS_seq);
		push(@length_ITS,$length_ITS);
		if ($length1 >= 100) {
			#if ($length1 < 1000) {
				if ($length2 >= 100) {
					#if ($length2 < 1000) {
						#print OUT2 ">$id\n$spacer1$ITS_seq$spacer2\n";
						@spacer1b = split(//,$spacer1);
						$length1b = scalar(@spacer1b);
						if ($length1b > 300) {
							$length1b_ = $length1b."***";
						}
						else {
							$length1b_ = $length1b;
						}
						push(@length1b,$length1b);
						@spacer2b = split(//,$spacer2);
						$length2b = scalar(@spacer2b);
						if ($length2b > 300) {
							$length2b_ = $length2b."***";
						}
						else {
							$length2b_ = $length2b;
						}
						push(@length2b,$length2b);
						@ITS_seqb = split(//,$ITS_seq);
						$length_ITS = scalar(@ITS_seqb);
						push(@length_ITSb,$length_ITS);
						print OUT2 ">$id|$length1b_|$length_ITS|$length2b_\n$spacer1$ITS_seq$spacer2\n";
					#}
				}
			#}
		}
		else {
			$i++;
			next;
		}
	}
	$i++;
}
close OUT2;

$min_length1 = min(@length1);
$min_length_ITS = min(@length_ITS);
$min_length2 = min(@length2);

$max_length1 = max(@length1);
$max_length_ITS = max(@length_ITS);
$max_length2 = max(@length2);

$mean_length1 = mean(@length1);
$mean_length_ITS = mean(@length_ITS);
$mean_length2 = mean(@length2);

open (OUT3,">>","ITS.stats") || die ("Error cannot write to ITS.stats: $!\n");
print OUT3 "ITS stats for all AFTOL sequences in Genbank with annotated spacer1, 5.8S, and spacer2 regions (before filtering for length)\n\n";
print OUT3 "Feature\tMin\tMax\tMean\n";
print OUT3 "ITS1:\t$min_length1\t$max_length1\t$mean_length1\n";
print OUT3 "5.8S:\t$min_length_ITS\t$max_length_ITS\t$mean_length_ITS\n";
print OUT3 "ITS2:\t$min_length2\t$max_length2\t$mean_length2\n";

$min_length1b = min(@length1b);
$min_length_ITSb = min(@length_ITSb);
$min_length2b = min(@length2b);

$max_length1b = max(@length1b);
$max_length_ITSb = max(@length_ITSb);
$max_length2b = max(@length2b);

$mean_length1b = mean(@length1b);
$mean_length_ITSb = mean(@length_ITSb);
$mean_length2b = mean(@length2b);

print OUT3"\n\nafter filtering for length\n\n";
print OUT3 "Feature\tMin\tMax\tMean\n";
print OUT3 "ITS1:\t$min_length1b\t$max_length1b\t$mean_length1b\n";
print OUT3 "5.8S:\t$min_length_ITSb\t$max_length_ITSb\t$mean_length_ITSb\n";
print OUT3 "ITS2:\t$min_length2b\t$max_length2b\t$mean_length2b\n";

close OUT3;
