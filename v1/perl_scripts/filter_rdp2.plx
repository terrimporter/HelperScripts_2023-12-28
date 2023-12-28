#!/usr/bin/perl
# Teresita M. Porter, Sept. 3, 2020
# Script to parse hmm output, grab Zotus, filter rdp.csv, then reformat for biom
# Usage perl filter_rdp2.plx hmm.txt rdp.out.tmp cutoff

use strict;
use warnings;
use Data::Dumper;
use Statistics::Descriptive;

# declare var
my $i=0;
my $line;
my $id;
my $score;
my $size;
my $stat;
my $percentile25;
my $percentile75;
my $iqr;
my $lc;
my $uc;
my $j;
my $seq;
my $record;
my $sk;
my $k;
my $p;
my $c;
my $o;
my $f;
my $g;
my $s;
my $sk_bp;
my $k_bp;
my $p_bp;
my $c_bp;
my $o_bp;
my $f_bp;
my $g_bp;
my $s_bp;
my $outfile = "obs_md.txt";
my $cutoff = $ARGV[2];
chomp $cutoff;
my $taxonomy;
my $confidence;

# declare array
my @hmm;
my @rdp;
my @line;
my @scores;
my @ids;

# declare hash
my %hmm; # key=id, value=score
my %keep; # key=id, value = 1

open (IN1, "<", $ARGV[0]) || die "Error can't open infile1: $!\n";
@hmm=<IN1>;
close IN1;

open (IN3, "<", $ARGV[1]) || die "Error can'open infile3: $!\n";
@rdp=<IN3>;
close IN3;

# parse HMMER output and hash scores
while ($hmm[$i]) {
	$line = $hmm[$i];
	chomp $line;

	if ($line =~ /^#/) {
		$i++;
		next;
	}
	else {
		@line = split ' ', $line; # split on whitespace (of any size)
		$id = $line[2];
		$score = $line[5];
		push(@scores, $score);
		$hmm{$id} = $score;

		$i++;
		next;
	}
}
$i=0;

# figure out cutoff for outlier hmmer scores
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@scores);
$percentile25 = $stat -> percentile(25);
$percentile75 = $stat -> percentile(75);
$iqr = $percentile75-$percentile25;
$lc = $percentile25-($iqr*1.5);
$uc = $percentile75+($iqr*1.5);

# get list of good ids (skip over ids with outlier scores)
while ( ($id, $score) = each(%hmm) ) {
	unless ($score < $lc || $score > $uc) {
		$keep{$id} = 1;
	}
}

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

# filter rdp by filtered ids
while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$id = $line[0];

	if (exists $keep{$id}) {
	
		parse_record();

		if ($i==0) { # headers
			print OUT "#OTUID\ttaxonomy\tconfidence\n";
		}

		print OUT $record."\t".$taxonomy."\t".$confidence."\n";

	}

	$i++;

}
$i=0;



#######################
sub parse_record {

	@line = split(/\t/, $line);

	$record = $line[0];

	$sk = $line[5];
	$k = $line[8];
	$p = $line[11];
	$c = $line[14];
	$o = $line[17];
	$f = $line[20];
	$g = $line[23];
	$s = $line[26];

	$sk_bp = $line[7];
	$k_bp = $line[10];
	$p_bp = $line[13];
	$c_bp = $line[16];
	$o_bp = $line[19];
	$f_bp = $line[22];
	$g_bp = $line[25];
	$s_bp = $line[28];

	if ($s_bp >= $cutoff){
		$confidence = $s_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k.";p__".$p.";c__".$c.";o__".$o.";f__".$f.";g__".$g.";s__".$s;
	}
	elsif ($g_bp >= $cutoff){
		$confidence = $g_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k.";p__".$p.";c__".$c.";o__".$o.";f__".$f.";g__".$g;
	}
	elsif ($f_bp >= $cutoff) {
		$confidence = $f_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k.";p__".$p.";c__".$c.";o__".$o.";f__".$f;
	}
	elsif ($o_bp >= $cutoff) {
		$confidence = $o_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k.";p__".$p.";c__".$c.";o__".$o;
	}
	elsif ($c_bp >= $cutoff) {
		$confidence = $c_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k.";p__".$p.";c__".$c;
	}
	elsif ($p_bp >= $cutoff) {
		$confidence = $p_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k.";p__".$p;
	}
	elsif ($k_bp >= $cutoff) {
		$confidence = $k_bp;
		$taxonomy = "root;sk__".$sk.";k__".$k;
	}
	elsif ($sk_bp >= $cutoff) {
		$confidence = $sk_bp;
		$taxonomy = "root;sk__".$sk;
	}
	else {
		$confidence = 0;
		$taxonomy = "Unassigned";
	}

}
