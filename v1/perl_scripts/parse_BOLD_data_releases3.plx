#!/usr/bin/perl
# Feb. 2/20 edited to get fasta and corresponding aa translation
#May 16, 2019 can't assume all BOLD seqs are Metazoa, need to map taxa higher than phyla to names and nodes.dmp
#July 4, 2018 edit to produce FASTA formatted file like testNBC.fasta
#April 30, 2018 parse all the BOLD data releases (tsv) into a FASTA file
#Automatically detect the different field names used to contain the same information
### Hardcoded cwd prefix
### Hardcoded list of files to read in, update these as needed
### Hardcoded, zero-indexed fields (column headers), update as needed
#USAGE perl parse_BOLD_data_releases2.plx 

use strict;
use warnings;

#declare variables
my $cwd="/home/terri/pseudogeneWork_2020/old2/";
my $outfile1="bold.nt.fasta";
my $outfile2 = "bold.aa.fasta";
my $namespath = "/home/terri/taxdmp_120419/names.dmp";
my $nodespath = "/home/terri/taxdmp_120419/nodes.dmp";
my $line;
my $i=0;
my $name;
my $nameclass;
my $taxid;
my $original_taxid;
my $taxid_list;
my $parent;
my $rank;
my $file;
my $fieldlist;
my $j=0;
my $header;
my $seq;
my $aaseq;
my $sampleID;
my $counter=0;
my $key;
my $value;
my $field;
my $field2;
my $n=0; #count of N's in seq (starting/trailing N's already removed, internal dashes changed to N's)
my $length;
my $aalength;
my $bin;
my $class;
my $phylum;
my $kingdom;
my $superkingdom;
my $count;

#declare array
my @names;
my @nodes;
my @line;
my @keys;
my @fieldlist;
my @in;
my @data;
my @value;
my @header;
my @seq;
my @aaseq;

#declare hash
my %taxid_name; #key=taxid, value=Genus species, values should be unique
my %name_taxid; #key=Genus species, value=taxid, can be a list of values
my %taxid_parent; #key=taxid, value=parent taxid
my %taxid_rank; #key=taxid, value=rank
my %headers; #key=data release file name, value=numeric fields to keep and the order to keep them
my %seqs; #key=data release file name, value=numeric field
my %aaseqs; #key=data release file name, value=numeric field
my %hash_header; #key=sampleID, value=header
my %hash_seq; #key=sampleID, value=seq
my %hash_aaseq; # key = sampleID, value = aaseq
my %hash_phyla; #key=phylum

open (NAMES, "<", $namespath) || die "Cannot open names.dmp: $!\n";
@names = <NAMES>;
close NAMES;

# hash names for quick lookups
while ($names[$i]) {
	$line = $names[$i];
	chomp $line;

	@line = split(/\t\|\t/,$line);
	$taxid = $line[0];
	$name = $line[1];
	$nameclass = $line[3];
	$nameclass =~ s/\t\|//g; #remove final pipe

	if ($nameclass eq "scientific name") {
		$taxid_name{$taxid} = $name; # should always be unique

		if (exists $name_taxid{$name}){ # name could have more than one taxid, do not overwrite
			$original_taxid = $name_taxid{$name};
			$taxid_list = join(";", $original_taxid, $taxid);
			$name_taxid{$name} = $taxid_list;
		}
		else {
			$name_taxid{$name} = $taxid;
		}
	}
	$i++;
}
$i=0;

open (NODES, "<", $nodespath) || die "Cannot open nodes.dmp: $!\n";
@nodes = <NODES>;
close NODES;

# hash nodes for quick lookups
while ($nodes[$i]) {
	$line = $nodes[$i];
	chomp $line;

	@line = split(/\t\|\t/,$line);
	$taxid = $line[0];
	$parent = $line[1];
	$rank = $line[2];

	$taxid_parent{$taxid} = $parent;
	$taxid_rank{$taxid} = $rank;

	$i++;
}
$i=0;

# identify which fields to keep for each data release version, zero-indexed fields
# field order: sampleid[1,2], accession[,35], phylum[8], class[9,10], order[10,12], family[11,14], genus[13,18], species[14,20], bin[,4]
#$headers{"CanadianBarcodeNet_ver1.txt"} = "2,,8,10,12,14,18,20,"; #inconsistent formatting noted in filename and file headers
$headers{"iBOL_phase_0.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_0.75_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_1.00_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_1.25_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_1.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_1.75_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_2.0_COI.tsv"}=       "1,33,8,9,10,11,13,14,4"; #inconsistent formatting noted in filename and file headers
$headers{"iBOL_phase_2.25_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_2.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_2.75_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase3.0_COI.tsv"}=        "1,34,8,9,10,11,13,14,4"; #inconsistent formatting noted in filename and file headers
$headers{"iBOL_phase_3.25_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_3.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_3.75_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_4.00_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_4.25_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_4.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4"; #problem with fields here to look at
$headers{"iBOL_phase_4.75_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_5.00_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_5.25_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_5.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_5.75_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_6.00_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_6.25_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";
$headers{"iBOL_phase_6.50_COI.tsv"}=      "1,35,8,9,10,11,13,14,4";

# identifiy which field contains the nucleotide sequence
#$seqs{"CanadianBarcodeNet_ver1.txt"} = "36"; #inconsistent formatting noted in filename and file headers
$seqs{"iBOL_phase_0.50_COI.tsv"}=      "30";
$seqs{"iBOL_phase_0.75_COI.tsv"}=      "30";
$seqs{"iBOL_phase_1.00_COI.tsv"}=      "30";
$seqs{"iBOL_phase_1.25_COI.tsv"}=      "30";
$seqs{"iBOL_phase_1.50_COI.tsv"}=      "30";
$seqs{"iBOL_phase_1.75_COI.tsv"}=      "30";
$seqs{"iBOL_phase_2.0_COI.tsv"}=       "30"; #inconsistent formatting noted in filename
$seqs{"iBOL_phase_2.25_COI.tsv"}=      "30";
$seqs{"iBOL_phase_2.50_COI.tsv"}=      "30";
$seqs{"iBOL_phase_2.75_COI.tsv"}=      "30";
$seqs{"iBOL_phase3.0_COI.tsv"}=        "30"; #inconsistent formatting noted in filename
$seqs{"iBOL_phase_3.25_COI.tsv"}=      "30";
$seqs{"iBOL_phase_3.50_COI.tsv"}=      "30";
$seqs{"iBOL_phase_3.75_COI.tsv"}=      "30";
$seqs{"iBOL_phase_4.00_COI.tsv"}=      "30";
$seqs{"iBOL_phase_4.25_COI.tsv"}=      "30";
$seqs{"iBOL_phase_4.50_COI.tsv"}=      "30";
$seqs{"iBOL_phase_4.75_COI.tsv"}=      "30";
$seqs{"iBOL_phase_5.00_COI.tsv"}=      "30";
$seqs{"iBOL_phase_5.25_COI.tsv"}=      "30";
$seqs{"iBOL_phase_5.50_COI.tsv"}=      "30";
$seqs{"iBOL_phase_5.75_COI.tsv"}=      "30";
$seqs{"iBOL_phase_6.00_COI.tsv"}=      "30";
$seqs{"iBOL_phase_6.25_COI.tsv"}=      "30";
$seqs{"iBOL_phase_6.50_COI.tsv"}=      "30";

# identifiy which field contains the amino acid sequence
#$seqs{"CanadianBarcodeNet_ver1.txt"} = ""; #inconsistent formatting noted in filename and file headers
$aaseqs{"iBOL_phase_0.50_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_0.75_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_1.00_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_1.25_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_1.50_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_1.75_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_2.0_COI.tsv"}=       "31"; #inconsistent formatting noted in filename
$aaseqs{"iBOL_phase_2.25_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_2.50_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_2.75_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase3.0_COI.tsv"}=        "31"; #inconsistent formatting noted in filename
$aaseqs{"iBOL_phase_3.25_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_3.50_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_3.75_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_4.00_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_4.25_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_4.50_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_4.75_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_5.00_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_5.25_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_5.50_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_5.75_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_6.00_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_6.25_COI.tsv"}=      "31";
$aaseqs{"iBOL_phase_6.50_COI.tsv"}=      "31";

#specify order to read in files, from oldest to newest
@keys = (
#"CanadianBarcodeNet_ver1.txt",
		"iBOL_phase_0.50_COI.tsv",
		"iBOL_phase_0.75_COI.tsv",
		"iBOL_phase_1.00_COI.tsv",
		"iBOL_phase_1.25_COI.tsv",
		"iBOL_phase_1.50_COI.tsv",
		"iBOL_phase_1.75_COI.tsv",
		"iBOL_phase_2.0_COI.tsv",
		"iBOL_phase_2.25_COI.tsv",
		"iBOL_phase_2.50_COI.tsv",
		"iBOL_phase_2.75_COI.tsv",
		"iBOL_phase3.0_COI.tsv",
		"iBOL_phase_3.25_COI.tsv",
		"iBOL_phase_3.50_COI.tsv",
		"iBOL_phase_3.75_COI.tsv",
		"iBOL_phase_4.00_COI.tsv",
		"iBOL_phase_4.25_COI.tsv",
		"iBOL_phase_4.50_COI.tsv",
		"iBOL_phase_5.00_COI.tsv",
		"iBOL_phase_5.25_COI.tsv",
		"iBOL_phase_5.50_COI.tsv",
		"iBOL_phase_5.75_COI.tsv",
		"iBOL_phase_6.00_COI.tsv",
		"iBOL_phase_6.25_COI.tsv",
		"iBOL_phase_6.50_COI.tsv");

#Add list of Eukaryote phyla from BOLD Taxonomy page
$hash_phyla{"Arthropoda"}="1";

#create an outfile
open (OUT, ">>", $outfile1) || die "Cannot open outfile: $!\n";
open (OUT2, ">>", $outfile2) || die "Cannot open outfile2: $!\n";

#open each file in order
while($keys[$i]) { #@keys contains filenames
	$file = $keys[$i];
	chomp $file;

	#identify the correct order of fields to grab from the data release files
	if (exists $headers{$file}) {
		$fieldlist = $headers{$file};
		@fieldlist = split(',',$fieldlist);

		#open file
		my $filename = $cwd.$file;
		open (IN, "<", $filename) || die "Cannot open infile $file: $!\n";
		@in = <IN>;

		#parse through file line by line
		while ($in[$j]) { #loop through data release file
			$line = $in[$j];
			chomp $line;

			#skip over header row
			if ($j == 0) {
				$j++;
				next;
			}
			elsif ($line =~ /WITHDRAWN/) {
				$j++;
				next;
			}
			else {
				#grab each field in data release file, automatically zero-indexed
				@data = split(/\t/,$line);

				#grab sampleID and phylum upfront to check for their presence in hashes
#				$field = $fieldlist[0];
#				$sampleID = $data[$field];
				$field2 = $fieldlist[2];
				$phylum = $data[$field2]; #start here because BOLD taxonomy doesn't go higher than this in the data releases!!!
				
				#if a value for phylum is found then proceed with processing
				if (length($phylum) > 0) {

					#subroutine to parse through each field and create a correctly formatted fasta header
					$header = parse_fields(); #pass to array references

					#if we have a Eukaryote phylum, continue processing, don't bother processing prokaryotes if present
					if (exists $hash_phyla{$phylum}) {
						$hash_header{$sampleID} = $header; #assign newly created header to sampleID (from sub)

						#verify which field contains the nucleotide sequence
						if (exists $seqs{$file}) {
							$field = $seqs{$file};
							$seq = $data[$field];

							# remove if sequence is missing
							if (! length($seq)) {
								$j++;
								next;
							}
							# remove starting N's or -'s
							if ($seq =~ m/^[N-]/) {
								$seq =~ s/^[N-]+//g;
								print "found leading N's\n";
							}
							# remove trailing N's or -'s
							if ($seq =~ m/[N-]$/) {
								$seq =~ s/[N-]+$//g;
								print "found trailing N's\n";
							}
							# change internal -'s to Ns
							if ($seq =~ m/-+/) {
								$seq =~ s/-/N/g;
								print "changed internal gaps to N's\n";
							}
							
							# dont't allow any ambiguous bases in seq, case insensitive search 
							unless ($seq =~ /[^ACGTacgt]/) { # try to match non-nucleotide bases
								@seq = split(//,$seq);
#								$length = scalar @seq;
								$hash_seq{$sampleID} = $seq;

								if (exists $aaseqs{$file}) {
									$field = $aaseqs{$file};
									$aaseq = $data[$field];

									if (length $aaseq) {
										@aaseq = split(//,$aaseq);
										$aalength = scalar @aaseq;
										$hash_aaseq{$sampleID} = $aaseq;
									}
									else {
										delete $hash_seq{$sampleID};
										$j++;
										next;
									}
								}
								else {
									delete $hash_seq{$sampleID};
									$j++;
									next;
								}
							}
							else {
								$j++;
								next;
							}	
						}
						else {
							print "can't find %seqs for $file\n";
						}
					}
					else {
#						print "phylum $phylum missing from hash_phylum\n";
					}
				}
			}
			$j++;
			$line=();
			@data=();
			$field=();
			$sampleID=();
			$class=();
			$phylum=();
			$kingdom=();
			$superkingdom=();
			$header=();
			$seq=();
			$aaseq=();
		}
		$j=0;
	}
	else {
		print "Cannot find file $file in headers hash\n";
	}
	$i++;
	$file=();
	$fieldlist=();
	@fieldlist=();
	@in=();
}
$i=0;

#print final non-redundant hash, all sampleIDs are unique
while (($key,$value) = each(%hash_header)) {
	$header = $value;
	$counter++;
	
	if (length $header) {
		if (exists $hash_seq{$key}) { #key=sampleID
			$seq = $hash_seq{$key};
			print OUT ">$header\n$seq\n";
		}
		if (exists $hash_aaseq{$key}) { #key = sampleID
			$aaseq = $hash_aaseq{$key};
			print OUT2 ">$header\n$aaseq\n";
		}
#		else {
#			print "sampleID $sampleID not found in hash_seq\n";	
#		}
	}
#	else {
#		print "Couldn't find seq for header $header\n";
#	}
}
#print "Counter=$counter\n";
close OUT;
close OUT2;

##############################################################################

sub parse_fields {

my $k=0;				
my $field;
my $value;
my $header;
my $newheader;
my $accession;

while ($fieldlist[$k]) {
	$field = $fieldlist[$k];
	chomp $field;
	if ($k==0) { # field forsampleID
		if (length $field) {
			$sampleID = $data[$field];
			chomp $sampleID; # need this in later loops
		}
		else {
			print "problem finding sampleID\n";
		}
	}
	elsif ($k==1) { # field for genbank accession
		if (length $field) { # field 1 and 8 could be empty
			$accession = $data[$field];
			if (length $accession) { #accession present
				if ($accession !~ m/^\d+/) { #make sure GB accession doesn't start with a number, probably a formatting error in the original file
					($kingdom,$superkingdom) = get_higher_level_fields();	
					$header = $accession." cellularOrganisms;".$superkingdom.";".$kingdom;
				}
				else {
					$accession="";
					$field = $fieldlist[8]; #BIN
					chomp $field;
					if (length $field) { # field 1 and 8 could be empty
						$bin = $data[$field];
						chomp $bin;
						if (length $bin) {	
							($kingdom,$superkingdom) = get_higher_level_fields();	
							$header = $bin." cellularOrganisms;".$superkingdom.";".$kingdom;
						}
						else {
							($kingdom,$superkingdom) = get_higher_level_fields();
							$header = $sampleID." cellularOrganisms;".$superkingdom.";".$kingdom;
						}
					}
					else {
						($kingdom,$superkingdom) = get_higher_level_fields();
						$header = $sampleID." cellularOrganisms;".$superkingdom.";".$kingdom;
					}
				}
			}
			else {
				$accession = "";
				$field = $fieldlist[8]; #BIN
				if (length $field) {
					$bin = $data[$field];
					chomp $bin;
					if (length $bin) {	
						($kingdom,$superkingdom) = get_higher_level_fields();	
						$header = $bin." cellularOrganisms;".$superkingdom.";".$kingdom;
					}
					else {
						($kingdom,$superkingdom) = get_higher_level_fields();
						$header = $sampleID." cellularOrganisms;".$superkingdom.";".$kingdom;
					}
				}
				else {
					($kingdom,$superkingdom) = get_higher_level_fields();				
					$header = $sampleID." cellularOrganisms;".$superkingdom.";".$kingdom;
				}
			}
		}
		else {
			$field = $fieldlist[8]; #BIN
			chomp $field;
			if (length $field) {
				$bin = $data[$field];
				chomp $bin;
				if (length $bin) {	
					($kingdom,$superkingdom) = get_higher_level_fields();	
					$header = $bin." cellularOrganisms;".$superkingdom.";".$kingdom;
				}
				else {
					($kingdom,$superkingdom) = get_higher_level_fields();
					$header = $sampleID." cellularOrganisms;".$superkingdom.";".$kingdom;
				}
			}
			else {
				($kingdom,$superkingdom) = get_higher_level_fields();				
				$header = $sampleID." cellularOrganisms;".$superkingdom.";".$kingdom;
			}
		}
	}
	elsif ($k > 1 && $k < 8) { #2 = phylum, #7 = species
		if (length $header) {
			$field = $fieldlist[$k];
			chomp $field;
			if (length $field) {
				$value = $data[$field];
				chomp $value;
				if (length $value) {
					$header = $header.";".$value;
				}
				else {
					$value="";
					$header = $header.";".$value;
				}
			}
			else {
				$value = "";
				$header = $header.";".$value;
			}
		}
		else {
			$value = "";
			$header = $header.";".$value;
		}
	}
	$k++;
}
return $header;
}
##############################################################################

sub get_higher_level_fields {

my $field;

	#for each phylum, check for kingdom and superkingdom
	$field = $fieldlist[2]; # phylum field
	$phylum = $data[$field]; 
	if (exists $name_taxid{$phylum}) {
		$taxid = $name_taxid{$phylum}; # need this for next sub
		$count = 0;
		($kingdom,$superkingdom) = check_level_up();
		return($kingdom,$superkingdom);
	}
	else {
		$phylum = "";
		$field = $fieldlist[3]; # class field
		$class = $data[$field];
		if ($class =~ /Oomycota/) {
			$class = "Oomycetes"; #BOLD uses Oomycota, NCBI uses Oomycetes at the class rank
		}
		if (exists $name_taxid{$class}) {
			$taxid = $name_taxid{$class}; #need this for next sub
			$count = 0;
			($kingdom,$superkingdom) = check_level_up();
			return($kingdom,$superkingdom);
		}
		else {
			$class = "";
			return("","");
		}
	}
}

##############################################################################

sub check_level_up {

	if (exists $taxid_parent{$taxid}) { # could be phylum or class rank taxid
		$parent = $taxid_parent{$taxid};
		if (exists $taxid_rank{$parent}) {
			$rank = $taxid_rank{$parent};
			if ($rank eq "kingdom") {
				$kingdom = $taxid_name{$parent};
				if ($count < 15) {
					$count++;
					$taxid = $parent;
					($kingdom,$superkingdom) = check_level_up();
				}
				else {
					return ($kingdom,"");
				}
			}
			elsif ($rank eq "superkingdom") {
				$superkingdom = $taxid_name{$parent};
				if (length $kingdom) {
					return($kingdom, $superkingdom);
				}
				else {
					$kingdom="";
					return($kingdom, $superkingdom);
				}
			}
			else {
				if ($count < 15) {
					$count++;
					$taxid = $parent;
					($kingdom,$superkingdom) = check_level_up();
				}
				else {
					return ("","");
				}
			}
		}
		else {
			return ("","");
		}
	}
	else {
		return ("","");
	}
}
