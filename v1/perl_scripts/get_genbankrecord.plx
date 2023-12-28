#!/usr/bin/perl
#January 24, 2010 by Terri Porter
#Updated Aug.31, 2010 to accept infile as an argument
#Use to retrieve genbank records using a list of gi's that were split into lists <=500 sequences using Split_ids_outfile.plx with outfiles called part1.out, part2.out, etc. that contains a vertical list of genbank gi numbers extracted from FastBlast_megablast.plx
#The outfile is called genbankrecords.outfile and needs to be renamed before running this script on another set of gi numbers
#USAGE: perl get_genbankrecord.plx part1.out

use Bio::DB::EUtilities; 
use Bio::SeqIO;   

#declare variables
my $line;
my $gi;
my $line1;
my $line2;
my $j=0; #flag accession
my $k=0; #flag id
my $lineage;

#declare arrays
my @ids;

#can use a list of gi numbers or gb numbers
open(IN,'<', $ARGV[0]) || die ("Error: $!\n");

while (<IN>) {
	$line = $_;
	chomp $line;
#	print $line."\n"; #test
	if ($line =~ /^\d+/){
		$line =~ /(\d+)/;
		$gi = $1;
		push (@ids, $gi);
	}
}
close(IN);
#print "@ids\n";#test

my $factory = Bio::DB::EUtilities->new(                          
	-eutil => 'efetch',                          
	-db => 'nucleotide',                          
	-rettype => 'genbank',                          
	-id => \@ids);   

my $file = 'genbankrecords.outfile';   
#dump HTTP::Response content to a file (not retained in memory) 
$factory->get_Response(-file => $file);	
