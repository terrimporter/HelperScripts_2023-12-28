#!/usr/bin/perl
#Dec. 9, 2016 script to parse RDP CO1v2 results by sample (year) for Blackbrook
#USAGE perl filter_rdp_by_site.plx 2014sites.txt 2015sites.txt BE_F230_concatenated.txt

#declare var
my $i=0;
my $line;
my $site;
my $header;

#declare array
my @sites2014; #blackbrook
my @sites2015; #blackbrook
my @rdp;
my @line;
my @header;

#declare hash
my %sites2014; #blackbrook
my %sites2015; #blackbrook

open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
@sites2014 = <IN>;
close IN;

foreach $site (@sites2014) {
	chomp $site;
	$sites2014{$site} = "1";
}

open (IN2, "<", $ARGV[1]) || die "Cannot open infile2: $!\n";
@sites2015 = <IN2>;
close IN2;

foreach $site (@sites2015) {
	chomp $site;
	$sites2015{$site} = "1";
}

open (IN3, "<", $ARGV[2]) || die "Cannot open infile3: $!\n";
@rdp = <IN3>;
close IN3;

open (OUT, ">>", "2014_rdp.out") || die "Cannot open outfile: $!\n";

open (OUT2, ">>", "2015_rdp.out") || die "Cannot open outfile2: $!\n";

while ($rdp[$i]) {
	$line = $rdp[$i];
	chomp $line;

	@line = split(/\t/,$line);

	$header = $line[0];
#	print "header: $header\n";#test
	@header = split(/;/,$header);
	$site = $header[1];
#	print "site: $site\n";#test

	if (exists $sites2014{$site}) {
		print OUT $line."\n";	
	}
	elsif (exists $sites2015{$site}) {
		print OUT2 $line."\n";
	}
	else {
		print "Couldn't assign OTU site $site to a file.\n";
	}
		
	$i++;
}
$i=0;
close OUT;
close OUT2;
