#!/usr/bin/perl
#April 26, 2010 by Terri Porter
#Script to parse sam file after using -N to disable iterative search with BWA, and using -n INT with samse when
#converting BWA .sai file to a .sam file, to count the number of transcripts mapped to
#USAGE perl parse_sam_XA.plx file.sam

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $flag;
my $transcriptID;
my $XAtag=();
my $j=0;
my $uniqueTranscripts;
my $transcriptIDline=();
my $k=0;

#declare array
my @in;
my @line;
my @XAtag;
my @transcriptIDline;

#declare hash
my %transcript;

open (IN,"<",$ARGV[0]) || die "Error cannot read sam file: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^@/) {
		$i++;
		next;
	}
	else {
		@line = split(/\t/,$line);
		#$readID = $line[0];
		$flag = $line[1];

		if ($flag != "4") { #only parse mapped reads
			#print "flag:$flag\n";#test
			$transcriptID = $line[2];
			#print "transcriptID:$transcriptID\n";#test
			$transcript{$transcriptID} = 1; #add to hash to easily remove duplicate keys/transcriptIDs
			
			$XAtag = $line[19];
			#print "XAtag:$XAtag\n";#test

			if ($XAtag) { #only parse XA tag if present
				$XAtag =~ s/XA:Z://;
				#print "XAtagline:$XAtag\n";
				@XAtag = split(/;/,$XAtag);

				while ($XAtag[$j]) {
					$transcriptIDline = $XAtag[$j];
					@transcriptIDline = split(/,/,$transcriptIDline);
					$transcriptID = $transcriptIDline[0];
					#print "transcriptID from line:$transcriptID\n";#test
			
					if (exists $transcript{$transcriptID}) {
						$j++;
						next;
					}
					else {
						$transcript{$transcriptID} = 1;
					}
					
					$j++;
					$transcriptIDline=();
					@transcriptIDline=();
					$transcriptID=();
				}
				$j=0;
				$XAtag=();
				@XAtag=();
			}
			$transcriptID=();
			$XAtag=();
		}
		@line=();
		$flag=();
	}
	$i++;
	$line=();
}
$i=0;

#test
#while(my($key,$value) = each (%transcript)) {
#	print "$key => $value\n";
#}

#grab unique transcriptIDs/keys and count them
$uniqueTranscripts = keys %transcript;

print "Reads were mapped to $uniqueTranscripts transcripts\n";
