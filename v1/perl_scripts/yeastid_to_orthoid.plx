#!/usr/bin/perl
#Aug. 30, 2012 by Terri Porter
#Script to use yeast.ids from grab_yeast_GOid.plx to grab orthoID from map.txt or SC_map.txt
#usage perl yeastid_to_orthoid.plx yeast.ids SC_map.txt

#declare var
my $i=0;
my $yeastID;
my $lengthYeastID;
my $zeroCount;
my $counter=0;
#my $j=0;
my $sgdid;
my $line;

#declare array
my @yeastID;
my @orthoID;
my @yeastID2;

#declare hash
my %sgdid;
my %orthoID;

open (YEASTID, "<", $ARGV[0]) || die "Error cannot open yeast.ids: $!\n";
@yeastID = <YEASTID>;
close YEASTID;

#print "@yeastID\n";#test

open (MAP, "<", $ARGV[1]) || die "Error cannot open map file: $!\n";
@orthoID = <MAP>;
close MAP;

#reformat yeastids into sgdid's
while ($yeastID[$i]) {
	$yeastID = $yeastID[$i];
	chomp $yeastID;

	@yeastID2 = split(//,$yeastID);
	$lengthYeastID = scalar(@yeastID2);
	$zeroCount = 9-$lengthYeastID;
	$sgdid = "SGDID:S";

	while ($counter < $zeroCount) {
		$sgdid = $sgdid."0";
		$counter++;
	}
	$counter=0;
	$sgdid = $sgdid.$yeastID;
	$sgdid{$sgdid} = 1;

	$i++;

	$yeastID=();
	@yeastID2=();
	$lengthYEastID=();
	$zeroCount=();
	$sgdid=();
}
$i=0;

#test
#while (my($key,$value) = each %sgdid) {
#	print $key."\n";
#}

#grab orthoIDs and sgdids
while ($orthoID[$i]) {
	$line = $orthoID[$i];
	chomp $line;

	if ($line =~ /^SC_\S+/) {
		$line =~ /^SC_(\S+)/;
		$orthoID = $1;
		if ($line =~ /SGDID:S\d+/) {
			$line =~ /(SGDID:S\d+)/;
			$sgdid = $1;
			$orthoID{$sgdid} = $orthoID;
		}
	}
	$i++;

	$line=();
	$orthoID=();
	$sgdid=();
}
$i=0;

open (OUT, ">>","yeast.orthoids") || die "Error cannot open outfile: $!\n";

while (($sgdid,$value) = each %sgdid) {
	$orthoID = $orthoID{$sgdid};
	print OUT $orthoID."\n";
	$orthoID=();
}
close OUT;
