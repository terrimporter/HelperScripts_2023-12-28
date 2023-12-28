#!/usr/bin/perl
# Teresita M. Porter, Sept. 2/21
# Script to add RDP classifier taxonomy to biom formatted ESV.biom (~ESV.table) from VSEARCH
# USAGE perl rdp_to_biom.plx rdp.out cutoff ESV.biom

use Data::Dumper;
use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $record; # ESV or OTU id
my $counter=1;
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
my $flag = 0;
my $cutoff = $ARGV[1];
chomp $cutoff;
my $confidence;
my $taxonomy;
my $outfile = "ESV_tax.biom";
my $shape;
my $id;
my $id_label;
my $new;
my $metadata;
my $metadata_label;
my $metadata_content;
my $j;
my $fix;

# declare array
my @rdp;
my @line;
my @id;
my @biom;
my @metadata;

# declare hash
my %confidence; # key = record id; value = confidence
my %taxonomy; # key = record id; value = QIIME-formatted taxonomy

# first parse rdp.out to filter each rank by cutoff
# second reformat using QIIME style
# third add to ESV.biom

open (IN, "<", $ARGV[0]) || die "Error cannot open rdp.out: $!\n";
@rdp = <IN>;
close IN;

# parse rdp assignments, filter by bootstrap support cutoff, hash
while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	parse_record();

	$confidence{$record} = $confidence;
	$taxonomy{$record} = $taxonomy;

	$i++;
	$confidence=();
	$taxonomy=();

}
$i=0;

open (IN2, "<", $ARGV[2]) || die "Error cannot open ESV.biom: $!\n";
@biom = <IN2>;
close IN2;

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

# process biom file
while ($biom[$i]) {
	$line = $biom[$i];
	chomp $line;

	if ($line =~ /\"shape\"/ && $flag == 0) {
		@line = split(',', $line);
		$shape = $line[0];
		$shape =~ s/\"shape\": \[//;
#		print "rows: ".$shape."\n"; #test
	}

	elsif ($line =~ /\"rows\"/ && $flag == 0) {
		$flag = 1;
		print OUT $line."\n";
		$i++;
		next;
	}

	elsif ($line =~ /\"id\"\:/ && $flag == 1 && $counter < $shape) {

		add_taxonomy();
		$counter++;

	}

	elsif ($line =~ /\"id\"\:/ && $flag == 1 && $counter == $shape) {

	add_taxonomy_to_last_record();	
	
	}

	elsif ($line =~ /\"columns\"/ && $flag == 1) {
		$flag=0;
		print OUT $line."\n";
		$i++;
		next;
	}

	else {
		print OUT $line."\n";
	}

	$i++;
	$confidence = ();
	$taxonomy = ();

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




#######################

sub add_taxonomy {

	@line = split(',', $line);
	$id = $line[0];
	@id = split(':', $id);
	$id_label = $id[0];
	$id_label =~ s/\{//;
	$record = $id[1];
	$record =~ s/\"//g;

	if (exists $confidence{$record}) {
		$confidence = $confidence{$record};
	}
	else {
		print "Can't find confidence for record $record\n";
	}

	if (exists $taxonomy{$record}) {
		$taxonomy = $taxonomy{$record};
	}
	else {
		print "Can't find taxonomy for record $record\n";
	}

	$new = "\n\t\t\t\t\"confidence\": \"$confidence\",\n\t\t\t\t\"taxonomy\": \"$taxonomy\"\n\t\t\t\}\n\t\t\},";

	$metadata = $line[1];
	@metadata = split(':', $metadata);
	$metadata_label = $metadata[0];
	$metadata_content = $metadata[1];
	$metadata_content =~ s/null}/$new/;
	$metadata = "\t\t\t".$metadata_label.": \{".$metadata_content;

	print OUT "\t\t\{\n\t".$id_label.":\"".$record."\",\n".$metadata."\n";

}




#######################

sub add_taxonomy_to_last_record {

	@line = split(',', $line);
	$id = $line[0];
	@id = split(':', $id);
	$id_label = $id[0];
	$id_label =~ s/\{//;
	$record = $id[1];
	$record =~ s/\"//g;

	if (exists $confidence{$record}) {
		$confidence = $confidence{$record};
	}
	else {
		print "Can't find confidence for record $record\n";
	}

	if (exists $taxonomy{$record}) {
		$taxonomy = $taxonomy{$record};
	}
	else {
		print "Can't find taxonomy for record $record\n";
	}

	$new = "\n\t\t\t\t\"confidence\": \"$confidence\",\n\t\t\t\t\"taxonomy\": \"$taxonomy\"\n\t\t\t\}\n\t\t\}";

	$metadata = $line[1];
	@metadata = split(':', $metadata);
	$metadata_label = $metadata[0];
	$metadata_content = $metadata[1];
	$metadata_content =~ s/null}/$new/;
	$metadata = "\t\t\t".$metadata_label.": \{".$metadata_content;

	print OUT "\t\t\{\n\t".$id_label.":\"".$record."\",\n".$metadata."\n";

}

