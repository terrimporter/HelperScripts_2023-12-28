#!/usr/bin/perl
#Jan. 6, 2012 by Terri Porter
#Script to make mothur .groups file, edited from make_mothur_infiles3.plx
#usage perl make_mothur_groups_file.plx seqtrim.fasta.noempty

#declare var
	my $i=0;
	my $seqtrim;
	my $line;
	my $id;
	my $group;
	my $groups_file;
	
	#declare array
	my @seqtrim;
	my @line;
	
	print "Please enter base name for .groups and .list files (ex. 16S):\n";
	$basename = <STDIN>;
	chomp $basename;
	$groups_file = $basename.".groups";
	
	open (OUT,">",$groups_file) || die ("Error cannot write to group file: $!\n");

#	while ($output[$i]) {
#		$seqtrim = $output[$i];
#		@seqtrim = split (/\./,$seqtrim);
#		$group = $seqtrim[0];
		
	#my $fastafile = $ARGV[0];
#		my @fastafile = split(/\./,$fastafile);
#		$group = $fastafile[0];

		open (IN,"<",$ARGV[0]) || die ("Error cannot read seqtrim file:$!\n");
		
		while (<IN>) {
			$line = $_;
			chomp $line;
			if (/>/) {
				@line = split(/ /,$line);
				$id = $line[0];
				$id =~ s/>//;
				print OUT "$id\t$basename\n";
			}
		}
#		close IN;
#		$i++;
#	}
	close OUT;

