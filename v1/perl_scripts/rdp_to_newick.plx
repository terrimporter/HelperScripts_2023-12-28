#!/usr/bin/perl
# Teresita M. Porter, June 8, 2019

# Script to turn RDP testNBC.taxonomy into Newick format to draw a tree in [r] ggtree read.tree
# Usage perl rdp_to_newick.plx testNBC.taxonomy

use strict;
use warnings;

# declare var
my $outfile = "taxonomy.nwk";


# edit these together (phylum tree use limit=3, switch=0; arthclass tree use limit=4, switch=1)
my $limit = 4; # rankid limit 3 - phylum, 4 - class
my $switch = 1; # 0 - all taxa, 1 - Arthropoda only

my $arthropoda = 0;
my $i=0;
my $line;
my $taxonid;
my $taxon;
my $parent;
my $rankid;
my $rank;
my $prevRank;

# declare array
my @in;
my @line;

# declare hash

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Error cannot open outfile: $!\n";

# testNBC.taxonomy
# 0-taxonid, 1-taxon, 2-parent taxid, 3-rankid, 4-rank name

# rankid
# 0-cellularOrganisms, 1-superkingdom, 2-kingdom, 3-phylum, 4-class, 5-order, 6-family, 7-genus, 8-species

# need to limit resolution of tree in the code itself

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\*/,$line);
	$taxonid = $line[0];
	$taxon = $line[1];
	$parent = $line[2];
	$rankid = $line[3];
	$rank = $line[4];
#	print $taxon."\t".$rank."\n"; #test

	if ($rankid == 0) { # first line in file
		print OUT "($taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}

	if ($rankid > $limit) { # limit the resolution of the tree so labels stay readable
		$i++;
		next;
	}

	if ($switch == 1) { # arthropoda only
		if ($rankid == 3) { # phylum
			if ($taxon eq "Arthropoda") {
				$arthropoda = 1;
				print_newick();
			}
			else {
				$arthropoda = 0;
				$i++;
				next;
			}
		}
		elsif ($rankid == 1) { # superkingdom
			print_newick();
		}
		elsif ($rankid == 2) { # kingdom
			print_newick();
		}
		else {
			if ($arthropoda == 1) {
				print_newick();
			}
			else {
				$i++;
				next;
			}
		}
	}
	else {
		print_newick();
	}

	$i++;
}
$i=0;

# add final closing brackets
while ($i < $prevRank+1) {
	print OUT ")";
	$i++;
}

# add final semicolon
print OUT ";\n";

close OUT;

#################################

sub print_newick {

	if ($rankid == $prevRank) { # same rank as last line
		print OUT ",$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}

	if ($rankid == $prevRank+1) { # finer rank
		print OUT ",($taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}

	if ($rankid == $prevRank-1) { # broader rank, 1 step
		print OUT ")$taxon";
		$prevRank = $rankid;
		$i++;
		next;

	}
	elsif ($rankid == $prevRank-2) { # broader rank, 2 step
		print OUT "))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	elsif ($rankid == $prevRank-3) { # broader rank, 3 step
		print OUT ")))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	elsif ($rankid == $prevRank-4) { # broader rank, 4 step
		print OUT "))))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	elsif ($rankid == $prevRank-5) { # broader rank, 5 step
		print OUT ")))))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	elsif ($rankid == $prevRank-6) { # broader rank, 6 step
		print OUT "))))))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	elsif ($rankid == $prevRank-7) { # broader rank, 7 step
		print OUT ")))))))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	elsif ($rankid == $prevRank-8) { # broader rank, 8 step
		print OUT "))))))))$taxon";
		$prevRank = $rankid;
		$i++;
		next;
	}
	else {
		print "Error, can't figure out how to process $line\n";
		# error here probably means no prevRank defined, ignore warning 
	}

}
