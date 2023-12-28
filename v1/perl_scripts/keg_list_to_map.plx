#!/usr/bin/perl
# Teresita M. Porter, Jan. 23/21
# Script to turn KEGG hierarchy table to a mapping file that is easier to work with
# KEGG mapfile available from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg
# USAGE perl keg_list_to_map.plx ko00001.keg

use strict;
use warnings;
use Data::Dumper;

# declare vars
my $outfile = "kegg.map";
my $i=0;
my $line;
my $id;
my $l1; # KEGG level 1
my $level;
my $l2; # KEGG level 2
my $l3; # KEGG level 3
my $l4; 

# declare arrays
my @in;
my @line;

# declare hashes

open (IN, "<", "ko00001.keg") || die "Error can't open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Error canpt open outfile: $!\n";
print OUT "KO\tLevel 1\tLevel 2\tLevel 3\tLevel 4\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^A/) { #level 1
		@line = split(/ /, $line);
		$id = shift @line;
		$l1 = join ' ', @line;
#		print "id: $id desc: $l1\n"; #test
	} 
	elsif ($line =~ /^B/) { # level 2
		@line = split(' ', $line);

		if (scalar @line > 1) {
			$level = shift @line;
			$id = shift @line;
			$l2 = join ' ', @line;
#			print "level: $level id: $id desc: $l2\n"; # test
		}
	}
	elsif ($line =~ /^C/) { # level 3
		@line = split(' ', $line);
		$level = shift @line;
		$id = shift @line;
		$l3 = join ' ', @line;
#		print "level: $level id: $id desc: $l3\n"; # test
	}
	elsif ($line =~ /^D/) { # level 4
		@line = split(' ', $line);
		$level = shift @line;
		$id = shift @line;
		$l4 = join ' ', @line;
#		print "level: $level id: $id desc: $l4\n"; #test
#		print "desc: $l4\n"; # test
		print OUT $id."\t".$l1."\t".$l2."\t".$l3."\t".$l4."\n";
	}

	$i++;

}
$i=0;
close OUT;


