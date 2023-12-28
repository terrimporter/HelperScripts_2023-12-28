#!/usr/bin/perl
#Feb. 6, 2017 updated how basename is retrieved
#Jan. 19, 2017 by Terri Porter
#Script to accomodate all the changes made to usearch v9.2 which no longer supports size= at every step and no longer enables easy clustering at values other than 97%; NOTE we use 98% for CO1 amplicons
#USAGE perl map_fasta_to_outtable.plx

use strict;
use warnings;

#declare variables
my $path_to_otus;
my $path_to_tables;
my $tablefile;
my $otufile;
my $path_to_otufile;
my $path_to_tablefile;
my $basename;
my $outfile;
my $i=0;
my $line;
my $otuID;
my $otusize;
my $NEWotusize;

#declare arrays
my @otu_files;
my @tables;
my @otufile;
my @tablefile;
my @line;
my @basename;

#declare hashes
my %tables; #key = basename, value = path to tablefile
my %otu_size; #key = otu, value = size

### Request that user enter directories containing files of OTUS and TABLES
print "Please enter path to the directory containing the filtered OTUs (including final '/'):\n";
$path_to_otus = <STDIN>;
chomp $path_to_otus;
#print $path_to_otus."\n"; #test

print "Please enter path to the directory containing the OTU tables (including final '/'):\n";
$path_to_tables = <STDIN>;
chomp $path_to_tables;
#print $path_to_tables."\n";#test

### Collect names for all OTU and TABLE files
opendir (OTUS, $path_to_otus) || die "Could not open $path_to_otus directory for reading: $!\n";
#skip files that are . or .. or other hidden files
#@otu_files = grep { /^\./ && -f "$path_to_otus/$_" } readdir(OTUS);
@otu_files = readdir(OTUS);
closedir OTUS;
#print @otu_files;#test

opendir (TABLES, $path_to_tables) || die "Could not open $path_to_tables directory for reading: $!\n";
#@tables = grep { /^\./ && -f "$path_to_tables/$_" } readdir(TABLES);
@tables = readdir(TABLES);
closedir TABLES;
#print @tables; #test

#hash @tables for easier searching
foreach $tablefile (@tables) {
	if ($tablefile =~ /^\./) {
		next;
	}
#	$basename = substr $tablefile,0,10;
	@basename = split(/\./,$tablefile);
	$basename = $basename[0];
	$path_to_tablefile = $path_to_tables.$tablefile;
	$tables{$basename} = $path_to_tablefile;
#	print "basename: $basename\n"; #test
#	print "path to tablefile: $path_to_tablefile\n"; #test
}

foreach $otufile (@otu_files) {
	chomp $otufile;

	if ($otufile =~ /^\./) {
		next;
	}

	$path_to_otufile = $path_to_otus.$otufile;
	open (IN, "<", $path_to_otufile) || die "Cannot open $otufile file: $!\n";
	@otufile = <IN>;
	close IN;

#	$basename = substr $otufile,0,10;
	@basename = split(/\./,$otufile);
	$basename = $basename[0];
	$outfile = $basename.".updatedsizefasta";

	open (OUT, ">>", $outfile) || die "Cannot open outfile $outfile: $!\n";

	if (exists $tables{$basename}) {
		$path_to_tablefile = $tables{$basename};

		open (IN2, "<", $path_to_tablefile) || die "Cannot open table $path_to_tablefile: $!\n";
		@tablefile = <IN2>;
		close IN2;

		#hash contents of table for easier lookups
		foreach $line (@tablefile) {
			chomp $line;
			
			if ($line =~ /^#/) {
				next;
			}
			else {
				@line = split (/\t/,$line);
				$otuID = $line[0];
				$otusize = $line[1];
				$otu_size{$otuID} = $otusize;
				@line = ();
				$otuID = ();
				$otusize = ();
			}
		}
			
		while ($otufile[$i]) {
			$line = $otufile[$i];
			chomp $line;

			if ($line =~ /^>/) {
				@line = split (/;/,$line);
				$otuID = $line[0];
				$otuID =~ s/^>//g;
				$otusize = $line[1];

				if (exists $otu_size{$otuID}) {
					$NEWotusize = $otu_size{$otuID};
					print OUT ">".$otuID."|size=".$NEWotusize."\n";
				}
				else {
					print "Error cannot find $otuID in table: $!\n";
				}

			}
			else {
				print OUT "$line\n";
			}
			$i++;
		}
		$i=0;
	}
	else {
		print "Error cannot find table for $basename:$!\n";
	}
	close OUT;
}

