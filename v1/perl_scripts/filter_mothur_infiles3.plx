#!/usr/bin/perl
#Aug.6,2013 edited to not print to STDOUT, too slow!
#Dec.22,2010 by Terri Porter
#Script to filter .groups and .list files from make_mothur_infiles.plx to remove singletons and doubletons prior to mothur analyses
#usage $perl filter_mothur_infiles.plx file.list file.groups

use strict;
use warnings;

#declare var
my $list_file;
my $line;
my $label;
my $num_otus;
my $i=0;
my $otu_line;
my $num_reads;
my $otu_line_to_keep;
my $x;
my $num_otus_filtered;
my $new_list_file;
my $new_otu_filtered_line;
my $group_file;
my $flag=0;
my $group_file_filtered;

#declare array
my @output;
my @line;
my @otu_line;
my @new_set;
my @reads_to_remove_from_group_file;
my @group_line_to_keep;
my @output2;


$list_file = $ARGV[0];

open (LIST,"<",$list_file) || die ("Error cannot open .list file: $!\n");

while (<LIST>) {
	$line = $_;
	chomp $line;
	@line = split(/\t/,$line);
	$label = shift(@line);
	$num_otus = shift(@line);
}
close LIST;

print "Closing $list_file.\n";
print "Removing singletons ...\n";

while ($line[$i]) {
	$otu_line = $line[$i];
	
	#print "Processing $otu_line\n";
	
	if ($otu_line =~ /,/) {
		@otu_line = split(/,/,$otu_line);
		$num_reads = scalar(@otu_line);
		if ($num_reads>=2) {###### modify read/otu cutoff here #####
			$otu_line_to_keep = join(',',@otu_line);
			push(@new_set,$otu_line_to_keep);
		}
		else {
			foreach $x (@otu_line) {
				push(@reads_to_remove_from_group_file, $x)
			}
		}
	}
	else {
		push(@reads_to_remove_from_group_file,$otu_line);
	}
	$i++;
}

$num_otus_filtered = scalar(@new_set);

print "There are $num_otus_filtered OTUs retained in the filtered set.\n";

$new_list_file = $list_file.".filtered";
open (NEWLIST,">>",$new_list_file) || die ("Error cannot write to new list file: $!\n");

$new_otu_filtered_line = join("\t",@new_set);

print NEWLIST "$label\t$num_otus_filtered\t$new_otu_filtered_line\n";
close NEWLIST;

@output2 = qx(ls | grep ".groups");
foreach $x (@output2) {
	print "$x\n";
}

print "Opening $group_file to process...\n";

$group_file = $ARGV[1];
open (GROUP,"<",$group_file) || die ("Error cannot open group file: $!\n");

while (<GROUP>) {
	$line = $_;
	chomp $line;
	
	#print "Processing $line\n";
	
	foreach $x (@reads_to_remove_from_group_file) {
		if ($line =~ /$x/) {
			$flag=1;
		}
	}
	if ($flag==0) {
		push(@group_line_to_keep,$line);
	}
	$flag=0;
}
close GROUP;

$group_file_filtered = $group_file.".filtered";
open (NEWGROUP,">>",$group_file_filtered) || die ("Error cannot write to new group file: $!\n");

foreach $x (@group_line_to_keep) {
	print NEWGROUP "$x\n";
}
close NEWGROUP;

print "New filtered .groups and .list files are ready for manual analysis in mothur.\n";
