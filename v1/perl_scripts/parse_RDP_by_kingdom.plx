#!/usr/bin/perl
#Dec.7, 2016 by Terri Porter
#Script to parse the output from RDP classifier CO1v2, filter by OTUsize and bootstrap support cutoff at certain rank
#USAGE perl parse_RDP_by_field.plx Athabasca_cat.txt

#declare var
my $cutoff;
my $minsize;
my $line;
my $family;
my $lengthsize;
my $size; #OTUsize, i.e. number of reads clustered in OTU
my $bp; #RDP classifier bootstrap proportion (0-1.0)

#declare array
my @parsed;
my @line;
my @lengthsize;

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@parsed = <IN>;
close IN;

print "Enter bootstrap support cutoff for family rank (0-1.0, CO1v2 use 0.3 for 200bp fragment):\n";
$cutoff = <STDIN>;
chomp $cutoff;

print "Enter minimum OTU cluster size cutoff (ex. minsize=3 will exclude singletons and doubletons and retain a minimum cluster size of 3 reads):\n";
$minsize = <STDIN>;
chomp $minsize;

print "Will use a bootstrap support cutoff of $cutoff and a minimum OTU size of $minsize\n";

open (OUT, ">>", "filtered_rdp.out") || die "Cannot open outfile: $!\n";

while ($parsed[$i]) {
	$line = $parsed[$i];
	chomp $line;
	@line = split(/\t/,$line);

	$family = $line[6]; #kingdom
#	print $family."\n";#check

	$lengthsize = $line[2];
#	print $lengthsize."\n";#check
	@lengthsize = split(/;/,$lengthsize);
	$size = $lengthsize[1];
	$size =~ s/size=//g;
#	print $size."\n";#check
	
	$bp = $line[7];
#	print $bp."\n";#check

	print "bp: $bp cutoff: $cutoff\n"; #test

	#check for minimum bootstrap support cutoff value
	if ($bp >= $cutoff) {
		#check for minimum OTU cluster size cutoff
		if ($size >= $minsize) {
			print OUT $line."\n";
		}
		else {
			$i++;
			next;
		}
	}
	else {
		$i++;
		next;
	}
	$i++;
}
$i=0;	
close OUT;
