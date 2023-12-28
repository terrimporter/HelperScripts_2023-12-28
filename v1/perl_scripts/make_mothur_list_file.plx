#!/usr/bin/perl
#Jan. 6, 2012 by Terri Porter
#Script to make mothur infiles, edited from make_mothur_infiles3.plx
#usage perl make_mothur_infiles4.plx results.clstr

#declare var
my $results_file;
my $label;
my $num_otus;
my $line;
my $id_part;
my $id;
my $group;
my $outfile1;
my $outfile2;
my $x;

#declare array
my @line;
my @id_part;
my @ids;

open (IN,"<",$ARGV[0]) || die "Error cannot read results.clstr file: $!\n";
my $clusterfile = $ARGV[0];

#####prompt for user input#####
print "Please enter label for this dataset (ex. uclust_097):\n";
$label = <STDIN>;
chomp $label;

$num_otus = qx(grep '^>' $clusterfile | wc -l);
chomp $num_otus;
print "num_otus: $num_otus\n";#test

$outfile1 = "seed.fasta.list";
open (OUT1,">>",$outfile1) || die "(Error cannot write to .list file: $!\n)";

print OUT1 "$label\t$num_otus";

while (<IN>) {
	$line = $_;
	chomp $line;
	if (/^>/) {
		if (/^>/) {
	#		@line = split (/\|/,$line);
	#		$id_part = $line[2];
#			@id_part = split(/ /,$id_part);
#			$id = $id_part[0];
#			push(@ids,$id);
#			print OUT1 "\t$id";
		}
	}
	else {
		if ($line =~ /\*/) {
			@line = split (/,/,$line);
			$id_part = $line[1];
			$id_part =~ s/\s+>//;
			$id_part =~ s/^>//;#add this
			$id_part =~ s/\.\.\.//;#add this
			@id_part = split(/ /,$id_part);
			$id = $id_part[0];
			#push(@ids,$id);
			#$line =~ />(\w+)\s+/;
			print OUT1 "\t$id";
		}
		else {
			@line = split (/,/,$line);
			$id_part = $line[1];
			$id_part =~ s/\s+>//;
			$id_part =~ s/^>//;#add this
			$id_part =~ s/\.\.\.//;#add this
			@id_part = split(/ /,$id_part);
			$id = $id_part[0];
			#push(@ids,$id);
			print OUT1 ",$id";
		}
	}
}

