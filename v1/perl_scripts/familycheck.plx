#!/usr/bin/perl
#Dec. 6, 2016 script to parse RDP CO1v2 results, parse out OTU assignment info for further processing
#USAGE perl familycheck.plx familylist.txt Athabasca_cat.txt

#declare var
my $i=0;
my $line;
my $site;
my $lengthsize;
my $size;
my $family;
my $familyBP;

#declare array
my @morphFamily; #families present in CABIN (morphology) set
my @rdp;
my @line;

#declare hash
my %morphFamily;

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@morphFamily = <IN>;
close IN;

foreach $morphFamily (@morphFamily) {
	chomp $morphFamily;
	$morphFamily{$morphFamily} = "1";
}

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@rdp = <IN2>;
close IN2;

open (OUT, ">>", "parsed.txt") || die "Cannot open outfile: $!\n";

while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	@line = split(/\t/,$line);

	$family = $line[14];
#	print $family,"\n";#test

	if (exists $morphFamily{$family}) {
		$site = $line[1];
		$lengthsize = $line[2];
	
		if ($lengthsize =~ m/size=\d+;/) {
			$lengthsize =~ m/size=(\d+)/;
			$size = $1;
		}
		else {
			print "Couldn't find OTU size in $site $lengthsize \n";
		}
	
		$familyBP = $line[15];

		print OUT $family."\t".$site."\t".$size."\t".$family."\t".$familyBP."\n";
		
	}
	else {
		$i++;
		next;
	}
	
	$i++;
}
$i=0;
