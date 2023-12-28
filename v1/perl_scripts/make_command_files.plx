#!/usr/bin/perl
#Sept. 11, 20122 edited to list all taxa as potential root
#Dec. 8, 2011 by Terri Porter
#Script to parse through RAxML_bootstrap.orhtoid files to get species order and potential outgroup
#also generate command.orthoid files to run phylip
#usage perl make_command_files.plx

use strict;
use warnings;

#declare var
my $dir;
my $i=0;
my $filename;
my $path_to_file;
my $j=0;
my $line;
my $taxon;
my $root;
my $flag=0;
my $root_id;

#declare array
my @dir;
my @file;
my @taxa;

#declare hash
my %taxa;

print "Enter directory path to bootstrap files including final / :\n";
$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Error cannot read from directory: $!\n";
@dir = readdir (DIR);
closedir DIR;

while ($dir[$i]) {
	$filename = $dir[$i];

	if ($filename =~ /\.\d+$/) {
		#print "file found\n";
			$path_to_file = $dir.$filename;
			
			open (IN, "<", $path_to_file) || die "Error cannot read from $filename: $!\n";
			@file = <IN>;
			close IN;
			
			$line = $file[0]; #get species order from first tree in file, just like phylip
			chomp $line;

			$line =~ s/\(//g;
			$line =~ s/\)//g; #remove brackets from newick file
			$line =~ s/;//;
			@taxa = split(/,/,$line); #keep in array to track order
			print "@taxa\n";#test

			while ($taxa[$j]) {  #put into hash to do quick searches
				$taxon = $taxa[$j];
				$taxa{$taxon} = $j;
				$j++;
			}
			$j=0;

			if ($flag==0 && $taxa{MB}) {
				$root = $taxa{MB}; #index in array
				$root_id = $root+1; #species number in list
				$flag=1;
				print "got MB root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{CA}) {
				$root = $taxa{CA};
				$root_id = $root+1;
				$flag=1;
				print "got AM root and passed conditional test\n";
				print_command_file();
			}

			elsif ($flag==0 && $taxa{AM}) {
				$root = $taxa{AM};
				$root_id = $root+1;
				$flag=1;
				print "got AM root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{BE}) {
				$root = $taxa{BE};
				$root_id = $root+1;
				$flag=1;
				print "got BE root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{BB}) {
				$root = $taxa{BB};
				$root_id = $root+1;
				$flag=1;
				print "got BB root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{GP}) {
				$root = $taxa{GP};
				$root_id = $root+1;
				$flag=1;
				print "got GP root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{PE}) {
				$root = $taxa{PE};
				$root_id = $root+1;
				$flag=1;
				print "got PE root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{NE}) {
				$root = $taxa{NE};
				$root_id = $root+1;
				$flag=1;
				print "got NE root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{BD}) {
				$root = $taxa{BD};
				$root_id = $root+1;
				$flag=1;
				print "got BD root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{SP}) {
				$root = $taxa{SP};
				$root_id = $root+1;
				$flag=1;
				print "got SP root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{MV}) {
				$root = $taxa{MV};
				$root_id = $root+1;
				$flag=1;
				print "got MV root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{EP}) {
				$root = $taxa{EP};
				$root_id = $root+1;
				$flag=1;
				print "got EP root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{PB}) {
				$root = $taxa{PB};
				$root_id = $root+1;
				$flag=1;
				print "got PB root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{MC}) {
				$root = $taxa{MC};
				$root_id = $root+1;
				$flag=1;
				print "got MC root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{RO}) {
				$root = $taxa{RO};
				$root_id = $root+1;
				$flag=1;
				print "got RO root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{SR}) {
				$root = $taxa{SR};
				$root_id = $root+1;
				$flag=1;
				print "got SR root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{UM}) {
				$root = $taxa{UM};
				$root_id = $root+1;
				$flag=1;
				print "got UM root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{LB}) {
				$root = $taxa{LB};
				$root_id = $root+1;
				$flag=1;
				print "got LB root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{SC}) {
				$root = $taxa{SC};
				$root_id = $root+1;
				$flag=1;
				print "got SC root and passed conditional test\n";
				print_command_file();
			}
			elsif ($flag==0 && $taxa{AN}) {
				$root = $taxa{AN};
				$root_id = $root+1;
				$flag=1;
				print "got AN root and passed conditional test\n";
				print_command_file();
			}


	}
	$i++;
	@file=();
	$line=();
	@taxa=();
	$taxon=();
	%taxa=();
	$flag=0;
	$root=();
}

sub print_command_file {

	#declare var
	my $orthoid;
	my $path_to_outfile;

	#declare array
	my @filename;

	@filename = split(/\./, $filename);
	$orthoid = $filename[1];
	$path_to_outfile = $dir.$orthoid.".command";
	open (OUT,">>", $path_to_outfile) || die "Error cannot write to command outfile: $!\n";
	print OUT "$filename\nC\nC\nC\nO\n$root_id\nY\n0.70\n";
	close OUT;

}

