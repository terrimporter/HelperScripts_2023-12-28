# Teresita M. Porter, Aug. 1, 2021
# Script to create a mapping file linking accession with lineage that was lost after VSEARCH processing
# watch out for missing ranks
# USAGE perl make_ncbi_mitofish_mapping.plx testNBC.fasta.ncbi testNBC.fasta.mitofish testNBC.taxonomy.mitofish

use strict;
use warnings;
use Data::Dumper;

# declare var
my $i=0;
my $line;
my $header;
my $acc;
my $lin; # lineage
my $num_ranks;
my $taxon;
my $prefix;
my $tax;
my $prev_taxon;
my $new_taxon;
my $new_lin;
my $outfile = "acc_tax.map";
my $rank;
my $r; # root
my $sk; # superkingdom
my $k;
my $p;
my $c;
my $o;
my $f;
my $g;
my $s; # species
my $lin2; # QIIME-formatted lineage

# declare array
my @ncbi;
my @header;
my @mftax; # MitoFish taxonomy file
my @mf; # MitoFish FASTA file
my @line;
my @lin; # lineage
my @taxon;
my @prefixes;

# declare hashes
my %map; # key = taxon, value = rank
my %rank = ("cellularOrganisms" => 'r__',
			"superkingdom" => 'sk__',
			"kingdom" => 'k__',
			"phylum" => 'p__',
			"class" => 'c__',
			"order" => 'o__',
			"family" => 'f__',
			"genus" => 'g__',
			"species" => 's__');
my %check; # key = prefix, value = prefix_taxon

open (IN, "<", $ARGV[0]) || die "Error can't open NCBI FASTA file: $!\n";
@ncbi = <IN>;
close IN;

open (OUT, ">>", $outfile) || die "Error can't open outfile: $!\n";

# parse ncbi FASTA headers in QIIME format
while ($ncbi[$i]) {
	$line = $ncbi[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/^>//g;

		@header = split(/\s+/, $header);
		$acc = $header[0];
		$lin = $header[1];
		@lin = split(/;/, $lin);
		$num_ranks = scalar(@lin);

		if ($num_ranks < 9) {
			check_ranks();
		}
		
		print OUT $acc."\t".$lin."\n";

	}
	$i++;
	%check=();
	$r=();
	$sk=();
	$k=();
	$p=();
	$c=();
	$o=();
	$f=();
	$g=();
	$s=();

}
$i=0;

# parse MitoFish FASTA headers
# edit to turn into QIIME format

open (IN2, "<", $ARGV[1]) || die "Error can't open MitoFish FASTA file: $!\n";
@mf = <IN2>;
close IN2;

open (IN3, "<", $ARGV[2]) || die "Error can't open MitoFish taxonomy file: $!\n";
@mftax = <IN3>;
close IN3;

# hash the taxonomy file for easier lookups
while ($mftax[$i]) {
	$line = $mftax[$i];
	chomp $line;

	@line = split(/\*/, $line);
	$taxon = $line[1];
	$rank = $line[4];
	$map{$taxon} = $rank;

	$i++;

}
$i=0;

#print Dumper(\%map);

while ($mf[$i]) {
	$line = $mf[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/^>//g;
		@header = split(/\t/, $header);
		$acc = $header[0];
		$lin = $header[1];

		@lin = split(/;/, $lin);

		# check each rank and compensate for any gaps
		check_ranks_reformat();
		
		$lin = join(";", $r, $sk, $k, $p, $c, $o, $f, $g, $s);
		print OUT $acc."\t".$lin."\n";

	}
	$i++;
	%check=();
	$r=();
	$sk=();
	$k=();
	$p=();
	$c=();
	$o=();
	$f=();
	$g=();
	$s=();

}
$i=0;

close OUT;


################

sub check_ranks {

	foreach $lin (@lin) {
		$taxon = $lin;
		@taxon = split(/__/, $taxon);
		$prefix = $taxon[0];
		$prefix = $prefix."__";
#		print "prefix $prefix\n";
		$tax = $taxon[1];
		$check{$prefix} = $taxon;
	}
		
	if (exists $check{"r__"}){
		$r = $check{"r__"};
	}
	else {
		print "Can't find root for $acc\n";
	}

	if (exists $check{"sk__"}) {
		$sk = $check{"sk__"};
	}
	else {
		$sk = "undef_".$r;
	}

	if (exists $check{"k__"}) {
		$k = $check{"k__"};
	}
	else {
		$k = "undef_".$sk;
	}

	if (exists $check{"p__"}) {
		$p = $check{"p__"};
	}
	else {
		$p = "undef_".$k;
	}

	if (exists $check{"c__"}) {
		$c = $check{"c__"};
	}
	else {
		$c = "undef_".$p;
	}

	if (exists $check{"o__"}) {
		$o = $check{"o__"};
	}
	else {
		$o = "undef_".$c;
	}

	if (exists $check{"f__"}) {
		$f = $check{"f__"};
	}
	else {
		$f = "undef_".$o;
	}

	if (exists $check{"g__"}) {
		$g = $check{"g__"};
	}
	else {
		$g = "undef_".$f;
	}

	if (exists $check{"s__"}) {
		$s = $check{"s__"};
	}
	else {
		$s = "undef_".$g;
		print "problem finding species for $acc\n";
	}

	$lin = join(";", $r, $sk, $k, $p, $c, $o, $f, $g, $s);

}



################
sub check_ranks_reformat {

	foreach $lin (@lin) {
		$taxon = $lin;
		if (exists $map{$taxon}) {
			$rank = $map{$taxon};
			if (exists $rank{$rank}) {
				$prefix = $rank{$rank};
				$taxon = $prefix.$taxon;
				$check{$rank} = $taxon;
			}
		}
	}

#	print Dumper(\%check);

	if (exists $check{"cellularOrganisms"}) {
		$r = $check{"cellularOrganisms"};
	}
	else {
		print "Can't find root for $acc\n";
	}

	if (exists $check{"superkingdom"}) {
		$sk = $check{"superkingdom"};	
	}
	else {
		$sk = "undef_".$r;
	}

	if (exists $check{"kingdom"}) {
		$k = $check{"kingdom"};	
	}
	else {
		$k = "undef_".$sk;
	}

	if (exists $check{"phylum"}) {
		$p = $check{"phylum"};	
	}
	else {
		$p = "undef_".$k;
	}

	if (exists $check{"class"}) {
		$c = $check{"class"};	
	}
	else {
		$c = "undef_".$p;
	}

	if (exists $check{"order"}) {
		$o = $check{"order"};	
	}
	else {
		$o = "undef_".$c;
	}

	if (exists $check{"family"}) {
		$f = $check{"family"};	
	}
	else {
		$f = "undef_".$o;
	}
		
	if (exists $check{"genus"}) {
		$g = $check{"genus"};	
	}
	else {
		$g = "undef_".$f;
	}

	if (exists $check{"species"}) {
		$s = $check{"species"};	
	}
	else {
		$s = "undef_".$g;
		print "problem finding species for $acc\n";
	}
		
	$lin = join(";", $r, $sk, $k, $p, $c, $o, $f, $g, $s);

}
