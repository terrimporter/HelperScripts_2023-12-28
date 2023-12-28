#!/usr/bin/perl	
#Terri Porter, June 3, 2013
#make list file for mothur using readid_OTUcentroid.map
#usage perl make_mothur_list.plx readid_OTUcentroid.map readid_sample.map

use strict;
use warnings;

#declare var
my $output;
my $i=0;
my $line;
my $readid;
my $OTUcentroid;
my $original;
my $new;
my $scalar;
my $scalar2;
my $max;
my $item;
my $flag=0; #flag
my $count=0;
my $sample;

#declare array
my @in;
my @line;
my @readid;
my @in2;

#declare hash
my %otu; #indexed by OTUcentroid
my %sample; #indexed by readid
my %otu_new;

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@in = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Can't open 2nd infile: $!\n";
@in2 = <IN2>;
close IN2;

open (OUT, ">>", "list.txt") || die "Cannot open outfile: $!\n";
print OUT "16Sv6_97otu\t";

#hash infile
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if (length($line) > 0) {
		@line = split(/\t/,$line);
		$readid = $line[0];
		$OTUcentroid = $line[1];
		if (exists $otu{$OTUcentroid}) {
			$original = $otu{$OTUcentroid};
			$new = $original."|".$readid;
			$otu{$OTUcentroid} = $new;
		}
		else {
			$otu{$OTUcentroid} = $OTUcentroid."|".$readid; #be sure to add OTUcentroid
		}
	}
	$i++;
	$line=();
	@line=();
	$readid=();
	$OTUcentroid=();
	$original=();
	$new=();
}
$i=0;

$scalar = keys(%otu);
print OUT "$scalar\t";

#hash 2nd infile
while($in2[$i]) {
	$line = $in2[$i];
	chomp $line;

	if (length($line) > 0 ) {
		@line = split(/\t/,$line);
		$readid = $line[0];
		$sample = $line[1];
		$sample =~ s/\D{1}$//g; #pool samples ABC, pool XZT soil cores
		if ($sample =~ /^F/) {
			$sample =~ s/^F_//g;
		}
		elsif ($sample =~ /^R/) {
			$sample =~ s/^R_//g;
		}
	}
	
	$sample{$readid} = $sample;
	$i++;
	$line=();
	@line=();
	$readid=();
	$sample=();
}
$i=0;

open (OUT2, ">>", "groups.txt") || die "Cannot open outfile2: $!\n";

#edit otu hash to remove reads that are not in the appropriate samples
while (($OTUcentroid,$readid) = each(%otu)) {
	@readid = split(/\|/,$readid);
	
	while ($readid[$i]) {
		$item = $readid[$i];
		$sample = $sample{$item};

		if ($sample =~ /(PAD3A|PAD3B|PAD3C)$/) {
			if (exists $otu_new{$OTUcentroid}) {
				$original = $otu_new{$OTUcentroid};
				$new = $original."|".$item;
				$otu_new{$OTUcentroid} = $new;
			}
			else {
				$otu_new{$OTUcentroid} = $item;
			}
		}
		$i++;
	}
	$i=0;
	@readid=();
}

while (($OTUcentroid,$readid) = each(%otu_new)) { ### add back in readid for the centroid itself!
	@readid = split(/\|/,$readid);
	$scalar2 = scalar(@readid);
	$max = $scalar2-1;

		while ($readid[$i]) {
			$item = $readid[$i];
#			$sample = $sample{$item};

			if ($i < $max) {
				$flag=1;
				print OUT "$item,";
				if (exists $sample{$item}) {
					$sample = $sample{$item};
					print OUT2 "$item\t$sample\n";
				}
				else {
					print "Can't find sample for $item\n";
				}
			}
			elsif ($i==$max) {
				$flag=1;
				print OUT "$item\t";
				if (exists $sample{$item}) {
					$sample = $sample{$item};
					print OUT2 "$item\t$sample\n";
				}
				else {
					print "Can't find sample for $item\n";
				}
			}

			$i++;
			$item=();
		}

		if ($flag == 1 ) {
			$count++;
		}

		$i=0;
		$flag=0;
	@readid=();
	$scalar2=();
	$max=();
}
close OUT;
close OUT2;

#print "Edit total number of OTUs by removing singletons $singletons\n";
print "Don't forget to change total number of OTUs to $count\n";
