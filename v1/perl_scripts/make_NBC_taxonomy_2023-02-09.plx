#!/usr/bin/perl
#Apr. 11, 2018 retain species rank for CO1 v3git training set; handle 'superkingdom' rank and create 'cellularOrganisms'
#Jan. 5, 2017 edit to retain species rank for CO1v3 training set
#Oct. 24, 2016 add superkingdom to hash
#Aug. 11, 2016 add Protura
#March 22, 2013 by Terri Porter
#Script to create taxonomy file for Ribosomal Database Project Naive Bayesian Classifier
#usage perl make_NBC_taxonomy.plx taxid.parsed.awk.uniq

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $cellularOrganisms;
my $superkingdom;
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
	
	$cellularOrganisms = $line[0];
	$cellularOrganisms =~ s/\s+//g; #removing trailing space
	
	$superkingdom = $line[1];
	$superkingdom =~ s/\s+//g; #remove any spaces
	
	$kingdom = $line[2];
	$kingdom =~ s/\s+//g;
	if ($kingdom eq 'undef') {
		$kingdom = $kingdom.'_'.$superkingdom;
	}
	$phylum = $line[3];
	$phylum =~ s/\s+//g;
	if ($phylum eq 'undef') {
		$phylum = $phylum.'_'.$kingdom;
	}
	
	$class = $line[4];
	$class =~ s/\s+//g;
	if ($class eq 'undef') {
		$class = $class.'_'.$phylum;
	}
	if ($class eq 'Deferribacteres') { # phylum and class of Bacteria
		$class = $class.'_'.$phylum;
	}

	$order = $line[5];
	$order =~ s/\s+//g;
	if ($order eq 'undef') {
		$order = $order.'_'.$class;
	}
	
	$family = $line[6];
	$family =~ s/\s+//g;
	if ($family eq 'undef') {
		$family = $family.'_'.$order;
	}
	
	$genus = $line[7];
	$genus =~ s/\s+//g;

	if ($genus eq 'undef') {
		$genus = $genus.'_'.$family;
	}
	
	$species = $line[8];
	$species =~ s/^\s+//g; #remove any preceeding whitespace
	$species =~ s/\s+/_/g; #change spaces to underscores just in case
	$species =~ s/_$//g; #if underscore at end of line, remove it
	if ($species eq 'undef') { # this should have already been removed with the make_NBC_fasta script
		$i++;
		next;
	}
	
	if (exists $lineage{$cellularOrganisms}) {
		
		if (exists $lineage{$superkingdom}) {
	
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
				create_kingdom();

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
		}
		else {
			create_superkingdom();

			if (exists $lineage{$kingdom}) {
				$i++;
				next;
			}
			else {
				create_kingdom();

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
		}
	}
	else {
		$termcounter++;
		$j=$termcounter-1;
		$lineage{$cellularOrganisms} = $termcounter."*".$cellularOrganisms."*".$j."*0*cellularOrganisms";
		
		if (exists $lineage{$superkingdom}) {
			$i++;
			next;
		}
		else {
			create_superkingdom();

			if (exists $lineage{$kingdom}) {
				$i++;
				next;
			}
			else {
				create_kingdom();

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
		}
	}
	$i++;
	$line=();
	@line=();
	$cellularOrganisms=();
	$superkingdom=();
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

open (OUT, ">>", "testNBC.taxonomy") || die "Error cannot open outfile:$!\n";

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
	$lineage{$species} = $termcounter."*".$species."*".$termcounter_previous."*8*species";

}

##########

sub create_genus {

	$termcounter++;
	$previous = $lineage{$family};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$genus} = $termcounter."*".$genus."*".$termcounter_previous."*7*genus";

}

###########

sub create_family {

	$termcounter++;
	$previous = $lineage{$order};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$family} = $termcounter."*".$family."*".$termcounter_previous."*6*family";

}

############

sub create_order {

	$termcounter++;
	$previous = $lineage{$class};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$order} = $termcounter."*".$order."*".$termcounter_previous."*5*order";

}

#############

sub create_class {

	$termcounter++;
	$previous = $lineage{$phylum};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$class} = $termcounter."*".$class."*".$termcounter_previous."*4*class";

}

##############

sub create_phylum {

	$termcounter++;
	$previous = $lineage{$kingdom};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$phylum} = $termcounter."*".$phylum."*".$termcounter_previous."*3*phylum";

}

###############

sub create_kingdom {

	$termcounter++;
	$previous = $lineage{$superkingdom};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$kingdom} = $termcounter."*".$kingdom."*".$termcounter_previous."*2*kingdom";

}

###############

sub create_superkingdom {

	$termcounter++;
	$previous = $lineage{$cellularOrganisms};
	@previous = split(/\*/, $previous);
	$termcounter_previous = $previous[0];
	$lineage{$superkingdom} = $termcounter."*".$superkingdom."*".$termcounter_previous."*1*superkingdom";

}
