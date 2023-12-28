#!/usr/bin/perl
#February 21, 2014 by Terri Porter
#Script to convert size in fasta header to percent.  Only work with 99% clusters.
#usage perl convert_size_to_percent.plx

use strict;
use warnings;
use List::Util 'sum';

#declare var
my $i=0;
my $file;
my $j=0;
my $line;
my $header;
my $id;
my $size;
my $sum;
my $well;
my $percent;
my $newfile;
my $flag=0;
my $length;

#declare array
my @output;
my @in;
my @line;

#declare hash
my %sum; #key = id, value = size
my %percent; #key = id, value = percent

@output = qx(ls | grep .centroids);

open (ID_SIZE, ">>", "id_size.map") || die "Error cannot open id_size.map: $!\n";

open (WELL_TOTALSIZE, ">>", "well_totalsize.map") || die "Error cannot open well_totalsize.map: $!\n";

open (ID_PERCENT, ">>", "id_percent.map") || die "Error cannot open id_percent.map: $!\n";

#parse through original fasta file with size in header
#create 3 map files:
# 1. id_size for all reads in each well
# 2. well_totalsize for all wells
# 3. id_percent for all reads in each well

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	open (IN, "<", $file) || die "Error, cannot open infile $file\n";
	@in = <IN>;
	close IN;

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;

		if ($line =~ /^>/) {
			$header = $line;
			
			if ($header =~ /^>\w+;size=\d+;/) {
				$header =~ /^>(\w+);size=(\d+);/;
				$id = $1;
				$size = $2;
				print ID_SIZE "$id\t$size\n";
				$sum{$id} = $size;

				$header=();
				$id=();
				$size=();
				
			}

			$header=();

		}

		$j++;
		$line=();

	}
	$j=0;

	$sum  = sum values %sum;
	$file =~ /^(\d+\w{1})/;
	$well = $1;
	print WELL_TOTALSIZE "$well\t$sum\n";

	while ( ($id,$size) = each (%sum) ) {
		$percent = sprintf '%.3f', $size / $sum;
		$percent{$id} = $percent;
		print ID_PERCENT "$id\t$percent\n";
		delete $sum{$id};
		$percent=();
	}

	%sum=();
	$sum=();
	$file=();
	$well=();
	@in=();

	$i++;

}
$i=0;

close ID_SIZE;
close WELL_TOTALSIZE;
close ID_PERCENT;

#parse through original fasta file, append with percent, but only print the ones with >= 0.250 percent.

while ($output[$i]) {
	$file = $output[$i];
	chomp $file;

	open (IN, "<", $file) || die "Error cannot open infile $file: $!\n";
	@in = <IN>;
	close IN;

	$newfile = $file.".percent";

	open (OUT, ">>", $newfile) || die "Error cannot open outfile $newfile: $!\n";

	while ($in[$j]) {
		$line = $in[$j];
		chomp $line;
#		$length = length ($line);

		if ($flag==0 && $line =~ /^>/) {
			$header = $line;

			if ($header =~ /^>\w+;size=\d+;/) {
				$header =~ /^>(\w+);size=\d+;/;
				$id = $1;

				if (exists $percent{$id}) {
					$percent = $percent{$id};

#					if ($percent >= 0) { ##### SS wants unfiltered OTUs with percent in fasta header (0.20) #####
					print OUT "$header";
					print OUT "percent=$percent;\n";
					$flag=1;
#					}
#					else {
#						$flag=0;
#					}
				}

			}
		}
		elsif ($flag==1 && $line !~ /^>/) {
			print OUT $line; #print all of seq on one line this time!
		}
		elsif ($flag==1 && $line =~ /^>/) {
			print OUT "\n"; 
			$flag=0; #reset
			$j--; #reprocess same line with $flag=0
		}

		$j++;
		$line=();
		$header=();
		$id=();
		$percent=();
	}
	$flag=0;#reset!!

	$j=0;

	$file=();
	$newfile=();
	close OUT;
	$i++;

}
$i=0;
