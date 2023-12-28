#!/usr/bin/perl
#Terri Porter, Sept. 8, 2016
#Cross-check target genera (or whatever taxon) against a concatenated report of all RDP taxonomic assignments
#USAGE perl xchecktargetgenera.plx targetgenera.sortedunique cat_rdp.txt

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $genus;
my $list;
my $id;
my $updatedlist;
my $target;
my $count;
my $bpp;

#declare array
my @targets;
my @assignments;
my @fields;
my @list;

#declare hash
my %genera; #keys are genera, values are pipe delimited SeqIDs

open (TARGET, "<", $ARGV[0]) || die "Cannot open first infile: $!\n";
@targets = <TARGET>;
close TARGET;

open (ASSTS, "<", $ARGV[1]) || die "Cannot open second infile: $!\n";
@assignments = <ASSTS>;
close ASSTS;

print "just read in two files\n";#test ok

my $test = scalar(@assignments);
print "array assignments: ".$test."\n";#test

#parse rdp taxonomic assignments into a hash with genus keys and pipe delimited seqid values
while ($assignments[$i]) {
	$line = $assignments[$i];
	chomp $line;
	@fields = split(/\t/, $line);
	$id = $fields[0]; #SeqID including sample, amplicon, OTUsize
	$genus = $fields[12]; #genus[12]
	$bpp = $fields[13]; #bootstrap proportion[13]
	$id = $id.$bpp.";"; #id now contains SeqID, sample, amplicon, OTUsize, GenusBP

	if (exists $genera{$genus}) {
		$list = $genera{$genus};
		$updatedlist = $list."|".$id;
		$genera{$genus} = $updatedlist;
	}
	else {
		$genera{$genus} = $id;
	}

	$i++;
	$line=();
	@fields=();
	$id=();
	$genus=();
	$list=();
	$updatedlist=();
}
$i=0;

open (OUT1, ">>", "target_genera_otu_count.txt") || die "Cannot open outfile1:$!\n";

open (OUT2, ">>", "target_genus_otuid.txt") || die "Cannot open outfile2:$!\n";

#for each target genus, check the hash, count OTUs, print a couple reports
while ($targets[$i]) {
	$target = $targets[$i];
	chomp $target;

	if (exists $genera{$target}) {
		$list = $genera{$target};
		
		if ($list =~ /\|/) {
			@list = split(/\|/, $list);	
			$count = scalar(@list);
			print OUT1 "$target\t$count\n";

			foreach $id (@list) {
				print OUT2 "$target\t$id\n";
			}

		}
		else {
			print OUT1 "$target\t1\n";
			print OUT2 "$target\t$list\n";#list here is just one value
		}
	}
	else {
		print OUT1 "$target\t0\n";
	}
	$i++;
	$target=();
	$list=();
	@list=();
	$count=();
	$id=();
}
$i=0;
close OUT1;
close OUT2;
