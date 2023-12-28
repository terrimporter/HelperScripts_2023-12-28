#!/usr/bin/perl
#Feb.15, 2013 edited to work with any fasta file, such as from Trimal that may have line breaks in seq part
#NEW USAGE perl generic_fasta_to_relaxed_phylip.plx file.fasta
#Nov.23,2010 by Terri Porter
#Script to convert fasta file (from mothur) into a relaxed phylip format for RAxML
#usage $perl fasta_to_relaxedphylip.plx mothur.fasta

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $line;
my $header;
my $i=0;
my $j;
my $nextline;
my $flag=0;
my $numseq=0;
my $char;
my $nextseq;
my $seq;
my $mean;
my $label;
my $scalar;
my $spaces;
#declare array
my @in;
my @seq;
my @char;
my @label;

#declare hash
my %label;
my %seq;

open (IN, "<", $ARGV[0]) || die "Error cannot open infile: $!\n";
@in = <IN>;
close IN;

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$header = $line;
		$header =~ s/>//;
		$numseq++;
		$label{$numseq} = $header;
		$header=();
	}
	elsif ($flag==1) { #need to concatenate seq from previous line
		$nextseq = $line;
		$seq = $seq.$nextseq;
		#check what's coming up on next line
		$j=$i+1;
		$nextline = $in[$j];
		chomp $nextline;

		if (defined $nextline && length $nextline > 0) {
			if ($nextline =~ /^>/) {
				$seq{$numseq} = $seq;
				$flag=0;
				@seq = split(//,$seq);
				$char = scalar (@seq);
				push(@char,$char);
			}
			elsif ($nextline !~ /^>/) {
				if ($nextline =~ /^\S+/) {
					$flag=1;	
				}
			}
		}
		else { #undefined, presumably end of file
			$seq{$numseq} = $seq;
		}
		$nextseq=();
	}
	else {
		$seq = $line;
		
		#check what's coming up on next line
		$j=$i+1;
		$nextline = $in[$j];

		if (defined $nextline && length $nextline > 0) {
			chomp $nextline;
			if ($nextline =~ /^>/) {
				$seq{$numseq} = $seq;
				$flag=0;
				@seq = split(//,$seq);
				$char = scalar (@seq);
				push(@char,$char);
			}
			elsif ($nextline !~ /^>/) {
				if ($nextline =~ /^\S+/) {
					$flag=1;	
				}
			}
		}
		else {
			print OUT "$seq\n";
			$seq{$numseq} = $seq;
			@seq = split(//,$seq);
			$char = scalar (@seq);
			push(@char,$char);
		}
	}
	$i++;
}
$i=0;

open (OUT, ">>", "relaxed.phy") || die "Error cannot open outfile: $!\n";

$mean = mean (@char);

print OUT "$numseq\t$mean\n";

foreach $numseq (sort keys %label) {
	$label = $label{$numseq};
	$label =~ s/\s{1}/_/g;
	@label = split(//,$label);
	$scalar = scalar(@label);
	$spaces = 30-$scalar;
	$seq = $seq{$numseq};
	print OUT "$label";
	print OUT " " x $spaces;
	print OUT "$seq\n";
	$spaces=();
	$scalar=();
}
close OUT;
