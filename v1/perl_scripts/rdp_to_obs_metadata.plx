#!/usr/bin/perl
# Teresita M. Porter, Sept. 3, 2021
# Convert RDP Classifier text output to a format that biom can use
# USAGE perl rdp_to_obs_metadata.plx rdp.out cutoff

use Data::Dumper;
use warnings;
use strict;

# declare vars
my $i=0;
my $line;
my $outfile = "obs_md.txt";
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

my $cutoff = $ARGV[1];

my $confidence;
my $taxonomy;

# declare arrays
my @in;
my @line;

# declare hashes


open (IN, "<", $ARGV[0]) || die "Error cannot open rdp.out: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	parse_record();

	if ($i==0) { # headers
		print OUT "#OTUID\ttaxonomy\tconfidence\n";
	}

	print OUT $record."\t".$taxonomy."\t".$confidence."\n";

	$i++;

}
$i=0;

close OUT;

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
#	print "s: $s\t"; #test

	$sk_bp = $line[7];
	$k_bp = $line[10];
	$p_bp = $line[13];
	$c_bp = $line[16];
	$o_bp = $line[19];
	$f_bp = $line[22];
	$g_bp = $line[25];
	$s_bp = $line[28];
#	print "s_bp: $s_bp\n"; #test

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
