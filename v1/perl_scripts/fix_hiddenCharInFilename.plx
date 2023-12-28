#!/usr/bin/perl
#
#



opendir (DIR, ".") ;
while ( my $file = readdir(DIR)) {
	chomp $file;
	if ($file =~ /1S\.merged$/) {
		print "|$file|\n";
		qx(mv $file 1S.merged_THISONE);
	}
#	if ($file eq "1S.merged\r\n") {
#		print "|$file|\n";
#	}
}

