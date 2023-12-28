#!/usr/bin/perl
# Teresita M. Porter, July 10, 2020
# Script to parse through FASTA files downloaded using BOLD API
# Filters for just COI, length 500-1500bp, no non-nucleotide characters, removes gaps
# Then looks up identification in names.dmp to grab taxid for next script to grab lineage
## 1. check if accession available, efetch gb, get taxid; if gb suppressed, gets taxid to higher rank, ok
## 2. check if identification in names.dmp, get taxid
## 3. check if first word in identification is genus in names.dmp, get taxid, create new dummy taxid to match identification (MARES ms Arranz et al., 2020)
## 4. if first word is not a proper genus, give up on it
### be sure to mention in ms that bold bin -> processid map would be helpful
# be sure to run dos2unix beforehand to change line endings for proper parsing!
# find . -type f -exec dos2unix {} \;
# USAGE perl filter_BOLD.plx

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

# vars

my $dir;
my $i=0;
my $filename;
my $j=0;
my $line;
my $processid;
my $identification;
my $marker;
my $accession;
my $k;
my $seq;
my $length;
my $outfile = "BOLD.fasta.filtered";
my $namedmp = "/home/terri/ncbi-blast-2.9.0+/db/names.dmp";
my $taxid;
my $name;
my $outfile2 = "BOLD.taxids";
my $temp;
my $words;
my $word;

# arrays
my @files;
my @fasta;
my @line;
my @seq;
my @names;
my @identification;

# hashes
my %names; #key = identification, value = taxid

print "Enter directory name containing BOLD.fasta files including final '/':\n";
$dir = <STDIN>;
chomp $dir;

opendir(DIR,$dir) || die "Cannot opendir $dir:$!\n";
@files = grep {/^[^\.]/ && -f "$dir/$_"} readdir(DIR);
close DIR;

open (IN, "<", $namedmp) || die "Cannot open names.dmp: $!\n";
@names = <IN>;
close IN;

# hash names.dmp for easier lookups
while ($names[$i]) {
	$line = $names[$i];
	chomp $line;

	@line = split(/\t\|\t/, $line);
	$taxid = $line[0];
	$name = $line[1];
	$names{$name} = $taxid;

	$i++;
}
$i=0;

open (my $fh_seq, ">>", $outfile) || die "Cannot open outfile: $!\n";
print $fh_seq "test\n";

open (my $fh_tax, ">>", $outfile2) || die "Cannot open outfile2: $!\n";
print $fh_tax "test2\n";

#=pod

while ($files[$i]) {
	$filename = $files[$i];
	chomp $filename;

	open (IN, "<", $dir."/".$filename) || die "Cannot open infile $filename: $!\n";
	@fasta = <IN>;
	close IN;

	while ($fasta[$j]) {
		$line = $fasta[$j];
		chomp $line;

		if ($line =~ /^>/) { #header
			@line = split(/\|/, $line);
			$processid = $line[0];
			$processid =~ s/^>//g;
			$identification = $line[1];
			$marker = $line[2];
			$accession = $line[3];

			if ($marker =~ /COI/) { #keep COI only
				$k = $j+1;
				$seq = $fasta[$k];
				chomp $seq;

				# get rid of any gaps 
				if ($seq =~ /\-/) {
					$seq =~ s/\-//g;
				}

				# if no non-nucleotide characters found then process
				if ($seq !~ /[^ACGT]/) {

					@seq = split(//, $seq);
					$length = scalar(@seq);

					# only process of length 500:1500 bp
					if ($length >= 500 && $length <= 1500) {

						if (defined $accession && length $accession > 0) { # accession is present

							if ($accession =~ /SUPPRESSED/) {
								$accession =~ s/\-SUPPRESSED//g;
							}
							#print $accession."\n";
							$taxid = get_taxid_from_gb_using_acc($accession);
							print "got taxid from gb file $taxid\n";
							print $fh_tax "$accession\t$taxid\n";
							print $fh_seq ">$accession\n$seq\n";

						}
						else {
							# no accesssion available
							# check if identification is in names.dmp

							if (exists $names{$identification}) {
								$taxid = $names{$identification};

								print $fh_tax "$processid\t$taxid\n";
								print $fh_seq ">$processid\n$seq\n";

							}
							else {
								# identification not in names.dmp
								# if multiple words, keep first (it may be a genus) and search names.dmp again

								@identification = split(/ /, $identification);
								$words = scalar(@identification);

								if ($words > 1) {
									$word = $identification[0];
									
									if (exists $names{$word}) {
										$taxid = $names{$word};

										print $fh_tax "$processid\t$taxid\n";
										print $fh_seq "$processid\n$seq\n";

									}
									else {
					#					print "Can't find word $word in names: $!\n";
									}
								}
								else {
										# only one word, give up on this.
								}
							}
						}
					}
					else {
						# seq too short or too long
					}
				}
				else {
					# non-nucleotide characters found
				}
			}
			else {
				# not COI
			}
		}
		else {
			# seq
		}
		$j++;
		$processid=();
		$identification=();
		$marker=();
		$accession=();
		$taxid=();
	}
	$j=0;
	$i++;
}
$i=0;

#=cut

close $fh_seq;
close $fh_tax;

#####

sub get_taxid_from_gb_using_acc {

	$accession = $_[0];

	# get gb record for accession
	my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
      	        	              	       -db      => 'nucleotide',
           		    	         		   -id      => $accession,
                   			      		   -email   => 'terriblue2002@yahoo.com',
                           			 	   -rettype => 'gbwithparts');

	$temp = $identification.".gb";
	print "temp $temp \n";

	# save response to file
	$factory -> get_Response(-file => $temp);

	# get taxon id from gb
	my $seqin = Bio::SeqIO -> new(	-file	=> $temp,
									-format	=> 'genbank');

	while (my $seq = $seqin -> next_seq) {
		my @feat_object = $seq -> get_SeqFeatures;

		foreach my $feat_object (@feat_object) {
			if ($feat_object -> primary_tag eq "source") {
				if ($feat_object -> has_tag('db_xref')) {
					my @value_dbxref = $feat_object -> get_tag_values('db_xref');
					foreach my $value_dbxref (@value_dbxref) {
						if ($value_dbxref =~ /taxon:/) {
							$value_dbxref =~ s/taxon://;
							$taxid = $value_dbxref;
								
							unlink $temp; #remove file when done
						}
					}
				}
			}
		}
	}
	$temp=();
	return $taxid;

}
