#!/usr/bin/perl
#Nov. 15, 2013 by Terri Porter
#Script to cross check taxa in tree file with taxa and samples in env.txt from create_unifrac_map.plx and category.txt from create_unifrac_cat.plx
#make sure raxml.tre is available in cwd
#usage perl crosscheck_unifrac_infiles_with_tree.plx env.txt category.txt

use strict;
use warnings;

#declare var
my $line;
my $readid;
my $sample;

#declare array
my @output;
my @env;
my @cat;
my @line;

#declare hash
my %readid; #index=readid, value = 1
my %sample; #index=sample, vlaue = 1

@output = system("tr , '\n' < raxml.tre > raxml.tre.edit");

print "Creating raxml.tre.edit\n";
open (TRE, "<", "raxml.tre.edit") || die "Error cannot open edited treefile: $!\n";

while (<TRE>) {
	$line = $_;
	chomp $line;

	if ($line =~ /\d+size\d+:/) {
		$line =~ /(\d+size\d+):/;
		$readid = $1;
#		print "readid:$readid\n";
		$readid{$readid} = 1;
	}
}

print "Creating env.txt.trimmed\n";
open (ENV, "<", $ARGV[0]) || die "Cannot open env file: $!\n";
open (ENVOUT, ">>", "env.txt.trimmed") || die "Cannot open env outfile: $!\n";

$readid=();
while (<ENV>) {
	$line = $_;
	chomp $line;

	@line = split(/\t/,$line);
	$readid = $line[0];
#	print "readid:$readid\n";
	$sample = $line[1];
	$sample{$sample} = 1;
#	print "sample:$sample\n";	
#	$abund = $line[2];

	if (exists $readid{$readid}) {
		print ENVOUT $line."\n";
	}

}
close ENV;
close ENVOUT;

print "Creating category.txt.trimmed\n";
open (CAT, "<", $ARGV[1]) || die "Cannot open category file: $!\n";
open (CATOUT, ">>", "category.txt.trimmed") || die "Error cannot open cat outfile: $!\n";

$sample=();
while (<CAT>) {
	$line = <CAT>;
	chomp $line;

	@line = split(/\t/,$line);
	$sample = $line[0];
#	$site = $line[1];

	if (exists $sample{$sample}) {
		print CATOUT $line."\n";
	}

}
close CAT;
close CATOUT;
