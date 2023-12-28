#!/usr/bin/perl
#Jan. 5, 2017 edit to retain species rank for CO1v3_species training set
#Aug. 11, 2016 add Protura
#March 22, 2013 by Terri Porter
#Script to create taxonomy file for Ribosomal Database Project Naive Bayesian Classifier
#usage perl make_NBC_taxonomy.plx taxid.parsed.species.awk.uniq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $kingdom;
my $phylum;
my $class;
my $order;
my $family;
my $genus;
my $species;
my $termcounter=0;
my $j;
my $key;
my $value;
my $previous;
my $termcounter_previous;

#declare array
my @in;
my @value;
my @line;
my @previous;

#declare hash
my %lineage;
my %sorted;

open (IN, "<", $ARGV[0]) || die "Error cannot open taxid.parsed.awk.uniq: $!\n";
@in = <IN>;
close IN;

#hash lineage table and create relationships
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	@line = split(/\t/, $line);
	$kingdom = $line[0];
	$kingdom =~ s/\s+//g; #remove any spaces
	$phylum = $line[1];
	$phylum =~ s/\s+//g;
	$class = $line[2];
	$class =~ s/\s+//g;
	$order = $line[3];
	$order =~ s/\s+//g;
	if ($order eq 'undef') {
		$order = $order.'_'.$class;
	}
	if ($order eq 'Plecoptera') { ### Plecoptera is genus of moths and order of stoneflies
		$order = $order.'_'.$class;
	}
	if ($order eq 'Protura') { ### Protura is class and order
		$order = $order.'_'.$class;
	}
	if ($order eq 'Diplura') { ### Diplura is class and order
		$order = $order.'_'.$class;
	}
	$family = $line[4];
	$family =~ s/\s+//g;
	if ($family eq 'undef') {
		$family = $family.'_'.$order;
	}
	$genus = $line[5];
	$genus =~ s/\s+//g;
	if ($genus eq 'undef') {
		$i++;
		next;
	}
	if ($genus eq 'Eutrapela') { ### Eutrapela is found in family Geometridae and Tenebrionidae
		$genus = $genus.'_'.$family;
	}
	if ($genus eq 'Plecoptera') {### Plecoptera is genus of moths and order of stoneflies
		$genus = $genus.'_'.$family;
	}
	$species = $line[6];
	$species =~ s/^\s+//g; #remove any trailing whitespace
	$species =~ s/\s+/_/g; #replace spaces with underscores just in case
	if ($species eq 'undef') {
		$i++;
		next;
	}

#	print "$kingdom\t$phylum\t$class\t$order\t$family\t$genus\n";

	if (exists $lineage{$kingdom}) {

		if (exists $lineage{$phylum}) {
		
			if (exists $lineage{$class}) {
			
				if (exists $lineage{$order}) {
					
					if (exists $lineage{$family}) {
						
						if (exists $lineage{$genus}) {

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
						else {
							create_genus();

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
					}
					else {
						create_family();

						if (exists $lineage{$genus}) {
							$i++;
							next;
						}
						else {
							create_genus();

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
					}
				}
				else {
					create_order();

					if (exists $lineage{$family}) {
						$i++;
						next;
					}
					else {
						create_family();

						if (exists $lineage{$genus}) {
							$i++;
							next;
						}
						else {
							create_genus();

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
					}
				}
			}
			else {
				create_class();

				if (exists $lineage{$order}) {
					$i++;
					next;
				}
				else {
					create_order();

					if (exists $lineage{$family}) {
						$i++;
						next;
					}
					else {
						create_family();

						if (exists $lineage{$genus}) {
							$i++;
							next;
						}
						else {
							create_genus();

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
					}
				}
			}
		}
		else {
			create_phylum();

			if (exists $lineage{$class}) {
				$i++;
				next;
			}
			else {
				create_class();

				if (exists $lineage{$order}) {
					$i++;
					next;
				}
				else {
					create_order();

					if (exists $lineage{$family}) {
						$i++;
						next;
					}
					else {
						create_family();

						if (exists $lineage{$genus}) {
							$i++;
							next;
						}
						else {
							create_genus();

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
					}
				}
			}
		}
	}
	else {
		$termcounter++;
		$j=$termcounter-1;
		$lineage{$kingdom} = $termcounter."*".$kingdom."*".$j."*0*kingdom";

		if (exists $lineage{$phylum}) {
			$i++;
			next;
		}
		else {
			create_phylum();

			if (exists $lineage{$class}) {
				$i++;
				next;
			}
			else {
				create_class();

				if (exists $lineage{$order}) {
					$i++;
					next;
				}
				else {
					create_order();

					if (exists $lineage{$family}) {
						$i++;
						next;
					}
					else {
						create_family();

						if (exists $lineage{$genus}) {
							$i++;
							next;
						}
						else {
							create_genus();

							if (exists $lineage{$species}) {
								$i++;
								next;
							}
							else {
								create_species();
							}
						}
					}
				}
			}
		}
	}
	$i++;
	$line=();
	@line=();
	$kingdom=();
	$phylum=();
	$class=();
	$order=();
	$family=();
	$genus=();
	$species=();
}
$i=0;
$termcounter=();

#sort the hash by termcounter into a new hash indexed by termcounter
while ( ($key, $value) = each (%lineage) ) {
	@value = split(/\*/, $value);
	$termcounter = $value[0];
	$sorted{$termcounter} = $value;
}
$key=();
$value=();

#sort keys then print

open (OUT, ">>", "testNBC_species.taxonomy") || die "Error cannot open outfile:$!\n";

foreach $key (sort {$a <=> $b} keys %sorted) {
	$value = $sorted{$key};
	print OUT $value."\n";
}
close OUT;

##########

sub create_species {

	$termcounter++;
	$previous = $lineage{$genus};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$species} = $termcounter."*".$species."*".$termcounter_previous."*6*species";

}

##########

sub create_genus {

	$termcounter++;
	$previous = $lineage{$family};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$genus} = $termcounter."*".$genus."*".$termcounter_previous."*5*genus";

}

###########

sub create_family {

	$termcounter++;
	$previous = $lineage{$order};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$family} = $termcounter."*".$family."*".$termcounter_previous."*4*family";

}

############

sub create_order {

	$termcounter++;
	$previous = $lineage{$class};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$order} = $termcounter."*".$order."*".$termcounter_previous."*3*order";

}

#############

sub create_class {

	$termcounter++;
	$previous = $lineage{$phylum};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$class} = $termcounter."*".$class."*".$termcounter_previous."*2*class";

}

##############

sub create_phylum {

	$termcounter++;
	$previous = $lineage{$kingdom};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$phylum} = $termcounter."*".$phylum."*".$termcounter_previous."*1*phylum";

}
