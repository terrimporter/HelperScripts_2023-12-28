#!/usr/bin/perl
#Nov. 14, 2013 by Terri Porter
#Create environment file for unifrac
#usage perl create_unifrac_map.plx readid_sample.map clusterreadid_readid.map otus.trimal.phy

use strict;
use warnings;

#declare var
my $line;
my $readid;
my $sample;
my $clusterreadid;
my $original;
my $new;
my $size;
my $readidlist;
my $value;
my $sizeline;

#declare array
my @line;
my @readidlist;
my @sizeline;

#declare hash
my %sample; #indexed = readid, value = sample
my %clusterreadid; #index = clusterreadid, value = readid
my %temp; #index = sample, value = readid++

open (MAP, "<", $ARGV[0]) || die "Error cannot open readid_sample.map: $!\n";

#hash readid-sample map file
print "Hashing readid-sample map...\n";
while(<MAP>) {
	$line = $_;
	chomp $line;
	@line = split(/\t/,$line);
	$readid = $line[0];
	$readid =~ s/ //g; #remove any spaces
	$sample = $line[1];
	$sample =~ s/ //g; #remove any spaces
	$sample{$readid} = $sample;

	$line=();
	@line=();
	$readid=();
	$sample=();

}
close MAP;

open (MAP2, "<", $ARGV[1]) || die "Error cannot open clusterreadid_readid.map: $!\n";

#hash clusterreadid-readid map file
print "Hashing clusterreadid-readid map...\n";
while(<MAP2>) {
	$line = $_;
	chomp $line;
	@line = split(/\t/,$line);
	$clusterreadid = $line[0];
	$clusterreadid =~ s/ //g;
	$readid = $line[1];
	$readid =~ s/ //g;

	if (exists $clusterreadid{$clusterreadid} ) {
		$original = $clusterreadid{$clusterreadid};
		$new = $original."|".$readid;
		$clusterreadid{$clusterreadid} = $new;
	}
	else {
		$clusterreadid{$clusterreadid} = $readid;
	}
	
	$line=();
	@line=();
	$clusterreadid=();
	$readid=();
	$original=();
	$new=();
}
close MAP2;

#create new environment file
open (OUT, ">>", "env.txt") || die "Error cannot open env.txt: $!\n";
open (IN, "<", $ARGV[2]) || die "Error cannot open otus.trimal.phy: $!\n";

print "Parsing phylip file to create infile for unifrac...\n";
my $i=1;
while (<IN>) {
	$line = $_;
	chomp $line;

	if ($line =~ /size/) {
#		print "line $i\n";
		@line = split(/_/,$line);
		$clusterreadid = $line[0];
		print "clusterreadid:$clusterreadid\n";
		$sizeline = $line[2];
		@sizeline = split(/ /,$sizeline);
		$size = $sizeline[0];

		if (exists $clusterreadid{$clusterreadid}) {
			$readidlist = $clusterreadid{$clusterreadid};
			@readidlist = split(/\|/,$readidlist);

			foreach $readid (@readidlist) {
				print "\treadid:$readid\n";
				if (exists $sample{$readid}) {
					$sample = $sample{$readid};
					print "\t\tsample:$sample\n";
					if (exists $temp{$sample}) {
						$original = $temp{$sample};
						$original++;
						$temp{$sample} = $original;
					}
					else {
						$temp{$sample} = 1; #initialize
					}
				}
				else {
					print "Cannot find sample for readid $readid\n";
				}
				$sample=();
				$original=();
			}

			while(($sample,$value) = each(%temp)) {
				print OUT $clusterreadid."_size_".$size."\t".$sample."\t".$value."\n";
			}
			%temp=();
		}
		else {
			print "Cannot find clusterreadid $clusterreadid\n";
		}
	}
	$line=();
	@line=();
	$clusterreadid=();
	$size=();
	$readidlist=();
	@readidlist=();
	$sizeline=();
	@sizeline=();
	$i++;
}
close IN;
close OUT;

